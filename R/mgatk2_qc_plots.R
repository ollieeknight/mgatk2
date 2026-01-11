library(tidyverse)
library(hdf5r)
library(Matrix)
library(Signac)
library(Seurat)

options(scipen = 999)

source("~/code/mgatk2/R/mgatk2_functions.R")

# Load data
mgatk_data <- read_mgatk_hdf5("code/mgatk2/tests/run_output/")

# Calculate cell coverage stats
cell_coverage_stats <- calculate_cell_coverage_stats(mgatk_data) %>%
  mutate(total_bases = mgatk_data$total_bases[barcode])

ggplot(cell_coverage_stats, aes(x = log10(mean_coverage), y = coverage_breadth)) +
  geom_point(alpha = 1, colour = 'darkblue') +
  labs(title = NULL,
       x = "Mean depth (log10)", 
       y = "Coverage breadth") +
  theme_classic() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = log10(10), linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.99, linetype = "dashed", color = "red")

# Filter cells
mgatk_data <- filter_cells_by_coverage(mgatk_data,
                                       min_mean_coverage = 10, 
                                       min_coverage_breadth = 0.99,
                                       cell_coverage_stats = cell_coverage_stats)

position_coverage <- tibble(
  position = mgatk_data$positions,
  mean_coverage = colMeans(mgatk_data$counts$coverage)
)

ggplot(position_coverage, aes(x = position, y = mean_coverage)) +
  geom_line(color = "darkblue", linewidth = 0.5) +
  labs(title = NULL, 
       x = "chrM (bp)", y = "Mean depth") +
  theme_classic() + 
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(mgatk_data$positions))) +
  scale_y_continuous(expand = c(0, 0))
  
ggplot(position_coverage, aes(x = position, y = mean_coverage)) +
  geom_line(color = "darkblue", linewidth = 0.5) +
  coord_polar() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(mgatk_data$positions))) +
  scale_y_continuous(expand = c(0, 0), trans = 'log10') +
  labs(title = NULL, 
       x = "chrM (bp)", 
       y = "Coverage") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.text.x = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    panel.grid.major = element_line(colour = alpha("black", 0.3)),
    panel.grid.minor = element_line(colour = alpha("black", 0.2))
  )

position_stats <- calculate_position_coverage_stats(mgatk_data)

ggplot(position_stats, aes(x = position, y = coverage_dropout)) +
  geom_line(color = "darkblue", linewidth = 0.8) +
  labs(title = NULL,
       x = "chrM (bp)", 
       y = "Fraction cells with 0 coverage") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(mgatk_data$positions))) +
  scale_y_continuous(expand = c(0, 0))
  
transposition_stats <- calculate_transposition_stats(mgatk_data)

ggplot(transposition_stats, aes(x = position, y = tn5_cuts_total)) +
  geom_col(fill = "darkblue", width = 10) +
  labs(title = NULL, 
       x = "chrM (bp)", y = "Tn5 cut sites (n)") +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

variants <- identify_variants(mgatk_data)

variants_filtered <- variants %>%
  filter(n_cells_conf_detected > 1)

ggplot(variants_filtered, aes(x = strand_correlation, y = vmr, 
                              colour = strand_correlation >= 0.65 & vmr > 0.1)) +
  geom_hline(yintercept = 0.10, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = 0.65, linetype = "dashed", colour = "black") +
  geom_point(alpha = 1) +
  scale_y_log10() +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "darkred")) +
  labs(
    title = NULL,
    x = "Strand concordance",
    y = "Variance-mean ratio"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1.01)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1))

variants_hq <- variants %>%
  filter(
    n_cells_conf_detected > 1,
    strand_correlation >= 0.65,
    vmr > 0.01
  )

allele_freq <- calculate_allele_freq(mgatk_data, variants_hq)

allele_freq_filtered <- allele_freq[apply(allele_freq, 1, max) >= 0.95, ]

# Visualize variant allele frequencies with pheatmap
library(pheatmap)

# Convert sparse matrix to regular matrix for pheatmap
allele_freq_mat <- as.matrix(allele_freq_filtered)

pheatmap(
  allele_freq_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = FALSE,
  color = colorRampPalette(c("white", "red"))(100),
  border_color = NA,
  fontsize_row = 8
  )

# seurat_object[['mgatk2']] <- CreateAssay5Object(allele_freq)

mito <- ReadMGATK(dir = 'code/mgatk2/tests/tenx_output/output/')

refallele <- mito$refallele

seurat_object <-  CreateSeuratObject(counts = mito$counts, assay = 'mito')

variable_sites <- IdentifyVariants(seurat_object, assay = "mito", refallele = refallele)
VariantPlot(variants = variable_sites, min.cells = 2)

confident_variants <- subset(variable_sites, 
                             subset = n_cells_conf_detected >= 2 & 
                               strand_correlation >= 0.65 & 
                               vmr > 0.01
)

seurat_object <- AlleleFreq(seurat_object, assay = 'mito', variants = confident_variants$variant)
seurat_object[['mgatk']] <- seurat_object[['alleles']]

DefaultAssay(seurat_object) <- 'mgatk'
seurat_object[['alleles']] <- NULL
seurat_object[['mito']] <- NULL

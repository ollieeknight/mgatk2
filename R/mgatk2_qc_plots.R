# Example QC and variant workflow for one mgatk2 run.
# From the repository root: Rscript R/mgatk2_qc_plots.R <mgatk2 output directory>

library(ggplot2)
source("R/mgatk2_functions.R")

dir <- commandArgs(trailingOnly = TRUE)[1]
data <- read_mgatk_hdf5(dir)

ggplot(data$cells, aes(log10(mean_depth), coverage_breadth)) +
  geom_point(colour = "darkblue") +
  labs(x = "Mean depth (log10)", y = "Coverage breadth") +
  theme_classic()

data <- filter_cells_by_coverage(data, min_mean_depth = 10, min_coverage_breadth = 0.9)
data <- recompute_reference_alleles(data)

positions <- position_stats(data)

ggplot(positions, aes(position, mean_depth)) +
  geom_line(colour = "darkblue") +
  labs(x = "chrM (bp)", y = "Mean depth") +
  theme_classic()

ggplot(positions, aes(position, tn5_fwd + tn5_rev)) +
  geom_col(fill = "darkblue", width = 10) +
  labs(x = "chrM (bp)", y = "Tn5 cut sites (n)") +
  theme_classic()

variants <- identify_variants(data, min_cells = 2)

ggplot(variants, aes(strand_correlation, vmr, colour = strand_correlation >= 0.65 & vmr > 0.01)) +
  geom_point() +
  geom_vline(xintercept = 0.65, linetype = "dashed") +
  geom_hline(yintercept = 0.01, linetype = "dashed") +
  scale_y_log10() +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "darkred"), guide = "none") +
  labs(x = "Strand correlation", y = "Variance-mean ratio") +
  theme_classic()

confident <- subset(variants, strand_correlation >= 0.65 & vmr > 0.01)
allele_freq <- calculate_allele_freq(data, confident)

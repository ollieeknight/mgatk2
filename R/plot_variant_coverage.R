library(tidyverse)
library(VariantAnnotation)

source("~/code/mgatk2/R/mgatk2_functions.R")

# Read VCF file
vcf <- read.table("~/code/mgatk2/data/mito_variants.vcf", 
                  comment.char = "#", 
                  header = FALSE,
                  col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))

variant_positions <- vcf$POS

# Read mgatk data
HC12_mgatk_data <- read_mgatk_hdf5('~/code/mgatk2/data/tapestri_HC12/')

# Calculate position coverage
HC12_position_coverage <- tibble(
  position = HC12_mgatk_data$positions,
  mean_coverage = colMeans(HC12_mgatk_data$counts$coverage),
  assay = "HC12",
  is_variant = position %in% variant_positions
)

# Plot 1: Full genome coverage with variants highlighted
ggplot(HC12_position_coverage, aes(x = position, y = mean_coverage)) +
  geom_line(color = "darkblue", linewidth = 0.5) +
  geom_point(data = HC12_position_coverage %>% filter(is_variant),
             aes(x = position, y = mean_coverage),
             color = "red", size = 3, shape = 21, fill = "red") +
  geom_vline(xintercept = variant_positions, 
             linetype = "dashed", color = "red", alpha = 0.3) +
  labs(title = "HC12 mitochondrial coverage with variants", 
       x = "chrM (bp)", 
       y = "Mean depth") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(HC12_mgatk_data$positions))) +
  scale_y_continuous(expand = c(0, 0), trans = 'log10')

# Plot 2: Zoomed in on variant positions only
variant_coverage <- HC12_position_coverage %>%
  filter(is_variant) %>%
  mutate(variant_label = paste0(position, "\nG>A"))

ggplot(variant_coverage, aes(x = factor(position), y = mean_coverage)) +
  geom_col(fill = "darkred", width = 0.7) +
  geom_text(aes(label = round(mean_coverage, 0)), 
            vjust = -0.5, size = 3.5) +
  labs(title = "Coverage at variant positions",
       x = "Position (chrM)",
       y = "Mean depth") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Plot 3: Coverage distribution across cells at each variant position
coverage_matrix <- HC12_mgatk_data$counts$coverage
variant_pos_indices <- which(HC12_mgatk_data$positions %in% variant_positions)

variant_coverage_per_cell <- map_dfr(variant_pos_indices, function(idx) {
  tibble(
    position = HC12_mgatk_data$positions[idx],
    coverage = coverage_matrix[, idx]
  )
})

ggplot(variant_coverage_per_cell, aes(x = factor(position), y = coverage)) +
  geom_violin(fill = "darkblue", alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  labs(title = "Coverage distribution at variant positions",
       x = "Position (chrM)",
       y = "Coverage per cell") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(trans = 'log10')

# Summary table
variant_summary <- HC12_position_coverage %>%
  filter(is_variant) %>%
  select(position, mean_coverage) %>%
  mutate(
    median_coverage = map_dbl(variant_pos_indices, ~median(coverage_matrix[, .x])),
    cells_covered = map_int(variant_pos_indices, ~sum(coverage_matrix[, .x] > 0)),
    total_cells = nrow(coverage_matrix)
  ) %>%
  arrange(position)

print(variant_summary)

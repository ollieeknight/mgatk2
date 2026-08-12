# Fixture check: build a tiny mgatk2 output, plant a variant, assert it comes back.
# Run with: Rscript R/test_mgatk2_functions.R

suppressPackageStartupMessages({
  library(hdf5r)
  library(dplyr)
  library(tibble)
})
source("R/mgatk2_functions.R")  # run from the repository root

n_positions <- 20
n_cells <- 12
variant_position <- 5
carriers <- 1:6

write_fixture <- function(dir) {
  dir.create(file.path(dir, "output"), recursive = TRUE)

  # hdf5r writes an R (cells x positions) matrix as HDF5 (positions x cells).
  empty <- function() matrix(0L, nrow = n_cells, ncol = n_positions)
  coverage <- matrix(20L, nrow = n_cells, ncol = n_positions)
  a_fwd <- matrix(10L, nrow = n_cells, ncol = n_positions)
  a_rev <- matrix(10L, nrow = n_cells, ncol = n_positions)
  g_fwd <- empty()
  g_rev <- empty()

  # Carrier i has 2i of its 20 reads on G, split evenly across the strands,
  # so heteroplasmy rises 0.1 to 0.6 and the strands correlate perfectly.
  a_fwd[carriers, variant_position] <- 10L - carriers
  a_rev[carriers, variant_position] <- 10L - carriers
  g_fwd[carriers, variant_position] <- carriers
  g_rev[carriers, variant_position] <- carriers

  counts <- H5File$new(file.path(dir, "output", "counts.h5"), mode = "w")
  counts[["barcode"]] <- sprintf("CELL%02d-1", seq_len(n_cells))
  counts[["A_fwd"]] <- a_fwd
  counts[["A_rev"]] <- a_rev
  counts[["G_fwd"]] <- g_fwd
  counts[["G_rev"]] <- g_rev
  for (name in c("C_fwd", "C_rev", "T_fwd", "T_rev", "tn5_cuts_fwd", "tn5_cuts_rev")) {
    counts[[name]] <- empty()
  }
  counts$close_all()

  reference <- rep("A", n_positions)
  reference[n_positions] <- "N"  # uncovered position, must yield no variants

  metadata <- H5File$new(file.path(dir, "output", "metadata.h5"), mode = "w")
  metadata[["coverage"]] <- coverage
  metadata[["reference"]] <- reference
  metadata[["mean_depth"]] <- rep(20, n_cells)
  metadata[["median_depth"]] <- rep(20, n_cells)
  metadata[["max_depth"]] <- rep(20L, n_cells)
  metadata[["genome_coverage"]] <- c(rep(1, n_cells - 2), 0.1, 0.1)
  metadata[["total_bases"]] <- rep(20 * n_positions, n_cells)
  metadata$close_all()
}

dir <- tempfile("mgatk2_fixture_")
write_fixture(dir)
data <- read_mgatk_hdf5(dir)

stopifnot(
  length(data$barcodes) == n_cells,
  length(data$positions) == n_positions,
  data$refallele[variant_position] == "A"
)

# Block reads must come back as cells x positions, in small blocks too.
block <- read_block(data, "A_fwd", 1:3)
stopifnot(nrow(block) == 3, ncol(block) == n_positions)
stopifnot(nrow(read_block(data, "A_fwd", 1L)) == 1)

variants <- identify_variants(data, min_cells = 1)
planted <- filter(variants, variant == paste0(variant_position, "A>G"))

heteroplasmy <- c(carriers / 10, rep(0, n_cells - length(carriers)))

stopifnot(
  nrow(planted) == 1,
  # Bulk frequency is alt reads over all reads, not the mean of the per-cell rates.
  isTRUE(all.equal(planted$mean, sum(2 * carriers) / (n_cells * 20))),
  planted$n_cells_detected == length(carriers),
  # Confident detection needs two reads on each strand, so carrier 1 misses out.
  planted$n_cells_conf_detected == sum(carriers >= 2),
  # Forward and reverse counts track each other exactly.
  isTRUE(all.equal(planted$strand_correlation, 1)),
  isTRUE(all.equal(planted$vmr, var(heteroplasmy) / (planted$mean + 1e-11))),
  # The reference base and the "N" position produce nothing.
  !any(variants$position == n_positions)
)

allele_freq <- calculate_allele_freq(data, planted)
stopifnot(
  dim(allele_freq) == c(1, n_cells),
  rownames(allele_freq) == paste0(variant_position, "A-G"),
  isTRUE(all.equal(unname(allele_freq[1, carriers]), carriers / 10)),
  all(allele_freq[1, -carriers] == 0)
)

position_stats <- calculate_position_coverage_stats(data)
stopifnot(
  nrow(position_stats) == n_positions,
  all(position_stats$mean_coverage == 20),
  all(position_stats$cells_10x == n_cells),
  all(position_stats$coverage_dropout == 0)
)

# Streaming across block boundaries must not change the answer.
blocked <- data
blocked$block_size <- 5
stopifnot(
  length(mgatk_blocks(blocked)) == 3,
  isTRUE(all.equal(variants, identify_variants(blocked, min_cells = 1))),
  isTRUE(all.equal(position_stats, calculate_position_coverage_stats(blocked)))
)

strand_stats <- calculate_strand_coverage_stats(data)
stopifnot(all.equal(strand_stats$mean_fwd_coverage, rep(10, n_positions)))

filtered <- filter_cells_by_coverage(data, min_mean_coverage = 10, min_coverage_breadth = 0.5)
stopifnot(length(filtered$barcodes) == n_cells - 2)

# Subsetting must survive a streaming read.
carrier_data <- subset_mgatk_barcodes(data, data$barcodes[carriers])
subset_freq <- calculate_allele_freq(carrier_data, planted)
stopifnot(isTRUE(all.equal(unname(subset_freq[1, ]), carriers / 10)))

close_mgatk_hdf5(data)
unlink(dir, recursive = TRUE)
cat("all checks passed\n")

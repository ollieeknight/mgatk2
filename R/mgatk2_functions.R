# Analysis helpers for mgatk2 HDF5 output.
# Matrices stay on disk; every function streams them in blocks of cells,
# which is the axis the files are chunked on.

read_mgatk_hdf5 <- function(mgatk_output_dir) {

  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required. Install with: install.packages('hdf5r')")
  }

  counts_file <- file.path(mgatk_output_dir, "output", "counts.h5")
  metadata_file <- file.path(mgatk_output_dir, "output", "metadata.h5")

  for (f in c(counts_file, metadata_file)) {
    if (!file.exists(f)) stop("File not found: ", f)
  }

  counts <- hdf5r::H5File$new(counts_file, mode = "r")
  metadata <- hdf5r::H5File$new(metadata_file, mode = "r")

  barcodes <- counts[["barcode"]][]
  per_cell <- function(name) setNames(metadata[[name]][], barcodes)

  barcode_metadata <- NULL
  if ("barcode_metadata" %in% names(metadata)) {
    group <- metadata[["barcode_metadata"]]
    fields <- setdiff(names(group), "barcode")
    barcode_metadata <- tibble(barcode = barcodes)
    for (field in fields) barcode_metadata[[field]] <- group[[field]][]
  }

  list(
    counts = counts,
    metadata = metadata,
    cells = seq_along(barcodes),
    barcodes = barcodes,
    positions = seq_len(metadata[["coverage"]]$dims[2]),
    refallele = metadata[["reference"]][],
    mean_depth = per_cell("mean_depth"),
    median_depth = per_cell("median_depth"),
    max_depth = per_cell("max_depth"),
    genome_coverage = per_cell("genome_coverage"),
    total_bases = per_cell("total_bases"),
    barcode_metadata = barcode_metadata
  )
}

close_mgatk_hdf5 <- function(mgatk_data) {
  mgatk_data$counts$close_all()
  mgatk_data$metadata$close_all()
  invisible(NULL)
}

# Blocks of cells. 512 keeps ~70 MB per matrix live; the files are chunked
# at 128 cells, so any multiple of 128 reads whole chunks. Override by setting
# mgatk_data$block_size.
mgatk_blocks <- function(mgatk_data) {
  size <- if (is.null(mgatk_data$block_size)) 512 else mgatk_data$block_size
  split(mgatk_data$cells, ceiling(seq_along(mgatk_data$cells) / size))
}

read_block <- function(mgatk_data, name, cells) {
  h5 <- if (name == "coverage") mgatk_data$metadata else mgatk_data$counts
  block <- h5[[name]][cells, ]
  if (is.null(dim(block))) matrix(block, nrow = length(cells)) else block
}

subset_mgatk_barcodes <- function(mgatk_data, barcodes) {

  keep <- which(mgatk_data$barcodes %in% barcodes)

  if (length(keep) == 0) stop("No matching barcodes found")

  message("  Found ", length(keep), " matching barcodes (",
          round(100 * length(keep) / length(mgatk_data$barcodes), 1), "%)")

  for (field in c("barcodes", "mean_depth", "median_depth", "max_depth",
                  "genome_coverage", "total_bases")) {
    mgatk_data[[field]] <- mgatk_data[[field]][keep]
  }
  if (!is.null(mgatk_data$barcode_metadata)) {
    mgatk_data$barcode_metadata <- mgatk_data$barcode_metadata[keep, ]
  }
  mgatk_data$cells <- mgatk_data$cells[keep]
  mgatk_data
}

filter_cells_by_coverage <- function(mgatk_data, min_mean_coverage, min_coverage_breadth) {

  keep <- mgatk_data$mean_depth >= min_mean_coverage &
          mgatk_data$genome_coverage >= min_coverage_breadth

  message("Kept ", sum(keep), " / ", length(keep), " cells")
  subset_mgatk_barcodes(mgatk_data, mgatk_data$barcodes[keep])
}

# mean/median/max depth and breadth are computed by the pipeline; only the
# depth-threshold columns need the coverage matrix.
calculate_cell_coverage_stats <- function(mgatk_data, deep = FALSE) {

  stats <- tibble(
    barcode = mgatk_data$barcodes,
    mean_coverage = unname(mgatk_data$mean_depth),
    median_coverage = unname(mgatk_data$median_depth),
    max_coverage = unname(mgatk_data$max_depth),
    coverage_breadth = unname(mgatk_data$genome_coverage)
  )

  if (!deep) return(stats)

  extra <- lapply(mgatk_blocks(mgatk_data), function(cells) {
    cov <- read_block(mgatk_data, "coverage", cells)
    tibble(
      coverage_cv = apply(cov, 1, sd) / rowMeans(cov),
      positions_10x = rowSums(cov >= 10),
      positions_50x = rowSums(cov >= 50)
    )
  })

  bind_cols(stats, bind_rows(extra))
}

# ponytail: no median_coverage here — a per-position median needs every cell
# resident at once. Use matrixStats::colMedians on a subset if you need it.
calculate_position_coverage_stats <- function(mgatk_data) {

  n_positions <- length(mgatk_data$positions)
  n_cells <- length(mgatk_data$cells)
  acc <- list(sum = numeric(n_positions), sum_sq = numeric(n_positions),
              max = numeric(n_positions), covered = numeric(n_positions),
              cells_10x = numeric(n_positions), cells_50x = numeric(n_positions))

  for (cells in mgatk_blocks(mgatk_data)) {
    cov <- read_block(mgatk_data, "coverage", cells)
    acc$sum <- acc$sum + colSums(cov)
    acc$sum_sq <- acc$sum_sq + colSums(cov^2)
    acc$max <- pmax(acc$max, apply(cov, 2, max))
    acc$covered <- acc$covered + colSums(cov > 0)
    acc$cells_10x <- acc$cells_10x + colSums(cov >= 10)
    acc$cells_50x <- acc$cells_50x + colSums(cov >= 50)
  }

  mean_coverage <- acc$sum / n_cells
  variance <- (acc$sum_sq - n_cells * mean_coverage^2) / (n_cells - 1)

  tibble(
    position = mgatk_data$positions,
    mean_coverage = mean_coverage,
    max_coverage = acc$max,
    coverage_cv = sqrt(pmax(variance, 0)) / mean_coverage,
    cells_covered = acc$covered,
    cells_10x = acc$cells_10x,
    cells_50x = acc$cells_50x,
    coverage_dropout = 1 - acc$covered / n_cells
  )
}

calculate_strand_coverage_stats <- function(mgatk_data) {

  n_positions <- length(mgatk_data$positions)
  fwd <- numeric(n_positions)
  rev <- numeric(n_positions)
  total <- numeric(n_positions)
  n_cells <- length(mgatk_data$cells)

  for (cells in mgatk_blocks(mgatk_data)) {
    for (base in c("A", "C", "G", "T")) {
      fwd <- fwd + colSums(read_block(mgatk_data, paste0(base, "_fwd"), cells))
      rev <- rev + colSums(read_block(mgatk_data, paste0(base, "_rev"), cells))
    }
    total <- total + colSums(read_block(mgatk_data, "coverage", cells))
  }

  tibble(
    position = mgatk_data$positions,
    mean_fwd_coverage = fwd / n_cells,
    mean_rev_coverage = rev / n_cells,
    total_coverage = total / n_cells
  ) %>%
    mutate(
      strand_balance = mean_fwd_coverage / pmax(total_coverage, 1),
      strand_imbalance = abs(mean_fwd_coverage - mean_rev_coverage)
    )
}

calculate_transposition_stats <- function(mgatk_data) {

  n_positions <- length(mgatk_data$positions)
  n_cells <- length(mgatk_data$cells)
  fwd <- numeric(n_positions)
  rev <- numeric(n_positions)
  cells_with_cuts <- numeric(n_positions)
  total_bases <- numeric(n_positions)

  for (cells in mgatk_blocks(mgatk_data)) {
    block_fwd <- read_block(mgatk_data, "tn5_cuts_fwd", cells)
    block_rev <- read_block(mgatk_data, "tn5_cuts_rev", cells)
    fwd <- fwd + colSums(block_fwd)
    rev <- rev + colSums(block_rev)
    cells_with_cuts <- cells_with_cuts + colSums(block_fwd + block_rev > 0)
    total_bases <- total_bases + colSums(read_block(mgatk_data, "coverage", cells))
  }

  tibble(
    position = mgatk_data$positions,
    tn5_cuts_total = fwd + rev,
    tn5_cuts_fwd = fwd,
    tn5_cuts_rev = rev,
    total_bases = total_bases,
    tn5_frequency = ifelse(total_bases > 0, (fwd + rev) / total_bases, 0),
    cells_with_tn5_cuts = cells_with_cuts,
    mean_tn5_cuts_per_cell = (fwd + rev) / n_cells,
    strand_bias = abs(fwd - rev) / pmax(fwd + rev, 1)
  )
}

# One streaming pass. Every statistic is a sum over cells, so nothing beyond
# per-position accumulators is ever resident.
identify_variants <- function(mgatk_data, min_cells = 0, min_strand_cor = 0, min_vmr = 0,
                              stabilise_variance = FALSE, low_coverage_threshold = 10) {

  bases <- c("A", "C", "G", "T")
  n_positions <- length(mgatk_data$positions)
  n_cells <- length(mgatk_data$cells)
  zeros <- function() matrix(0, nrow = n_positions, ncol = 4, dimnames = list(NULL, bases))
  acc <- list(alt = zeros(), af = zeros(), af_sq = zeros(), af_high = zeros(),
              af_sq_high = zeros(), n_low = zeros(), detected = zeros(), conf = zeros(),
              n = zeros(), f = zeros(), r = zeros(), f_sq = zeros(), r_sq = zeros(),
              fr = zeros())
  coverage_total <- numeric(n_positions)

  for (cells in mgatk_blocks(mgatk_data)) {
    cov <- read_block(mgatk_data, "coverage", cells)
    coverage_total <- coverage_total + colSums(cov)
    low <- cov < low_coverage_threshold

    for (base in bases) {
      fwd <- read_block(mgatk_data, paste0(base, "_fwd"), cells)
      rev <- read_block(mgatk_data, paste0(base, "_rev"), cells)
      alt <- fwd + rev
      af <- ifelse(cov > 0, alt / cov, 0)
      # Upstream takes cells carrying the alt on either strand, not just forward.
      seen <- fwd > 0 | rev > 0

      acc$alt[, base] <- acc$alt[, base] + colSums(alt)
      acc$af[, base] <- acc$af[, base] + colSums(af)
      acc$af_sq[, base] <- acc$af_sq[, base] + colSums(af^2)
      acc$af_high[, base] <- acc$af_high[, base] + colSums(af * !low)
      acc$af_sq_high[, base] <- acc$af_sq_high[, base] + colSums(af^2 * !low)
      acc$n_low[, base] <- acc$n_low[, base] + colSums(low)
      acc$detected[, base] <- acc$detected[, base] + colSums(alt > 0)
      acc$conf[, base] <- acc$conf[, base] + colSums(fwd >= 2 & rev >= 2)
      acc$n[, base] <- acc$n[, base] + colSums(seen)
      acc$f[, base] <- acc$f[, base] + colSums(fwd * seen)
      acc$r[, base] <- acc$r[, base] + colSums(rev * seen)
      acc$f_sq[, base] <- acc$f_sq[, base] + colSums(fwd^2 * seen)
      acc$r_sq[, base] <- acc$r_sq[, base] + colSums(rev^2 * seen)
      acc$fr[, base] <- acc$fr[, base] + colSums(fwd * rev * seen)
    }
  }

  results <- lapply(bases, function(base) {
    # "N" marks positions the pipeline never saw covered.
    idx <- which(mgatk_data$refallele != base & mgatk_data$refallele != "N")
    if (length(idx) == 0) return(NULL)

    bulk <- acc$alt[idx, base] / coverage_total[idx]

    if (stabilise_variance) {
      # Low-coverage cells are held at the bulk frequency, so their
      # contribution to both moments is known without a second pass.
      n_low <- acc$n_low[idx, base]
      sum_af <- acc$af_high[idx, base] + n_low * bulk
      sum_af_sq <- acc$af_sq_high[idx, base] + n_low * bulk^2
    } else {
      sum_af <- acc$af[idx, base]
      sum_af_sq <- acc$af_sq[idx, base]
    }
    variance <- (sum_af_sq - sum_af^2 / n_cells) / (n_cells - 1)

    n <- acc$n[idx, base]
    f <- acc$f[idx, base]
    r <- acc$r[idx, base]
    covariance <- n * acc$fr[idx, base] - f * r
    spread <- (n * acc$f_sq[idx, base] - f^2) * (n * acc$r_sq[idx, base] - r^2)
    strand_correlation <- ifelse(n > 1 & spread > 0, covariance / sqrt(spread), 0)

    tibble(
      position = mgatk_data$positions[idx],
      nucleotide = paste0(mgatk_data$refallele[idx], ">", base),
      variant = paste0(mgatk_data$positions[idx], mgatk_data$refallele[idx], ">", base),
      mean = ifelse(is.finite(bulk), bulk, 0),
      # Pseudocount matches Signac; bulk is 0 wherever nothing was observed.
      vmr = pmax(variance, 0) / (bulk + 1e-11),
      n_cells_detected = acc$detected[idx, base],
      n_cells_conf_detected = acc$conf[idx, base],
      strand_correlation = strand_correlation
    )
  })

  bind_rows(results) %>%
    filter(is.finite(vmr),
           n_cells_conf_detected >= min_cells,
           strand_correlation >= min_strand_cor,
           vmr > min_vmr) %>%
    arrange(desc(vmr))
}

# Dense variants x cells. Sparsify with as(x, "sparseMatrix") if you need it.
calculate_allele_freq <- function(mgatk_data, variants, min_coverage = 1) {

  allele_freq <- matrix(0, nrow = nrow(variants), ncol = length(mgatk_data$cells),
                        dimnames = list(gsub(">", "-", variants$variant), mgatk_data$barcodes))
  if (nrow(variants) == 0) return(allele_freq)

  alt_base <- substr(variants$nucleotide, 3, 3)
  offset <- 0

  for (cells in mgatk_blocks(mgatk_data)) {
    cov <- read_block(mgatk_data, "coverage", cells)[, variants$position, drop = FALSE]
    columns <- seq_along(cells) + offset

    for (base in unique(alt_base)) {
      rows <- which(alt_base == base)
      positions <- variants$position[rows]
      alt <- read_block(mgatk_data, paste0(base, "_fwd"), cells)[, positions, drop = FALSE] +
             read_block(mgatk_data, paste0(base, "_rev"), cells)[, positions, drop = FALSE]
      block_cov <- cov[, rows, drop = FALSE]
      allele_freq[rows, columns] <- t(ifelse(block_cov >= min_coverage, alt / block_cov, 0))
    }
    offset <- offset + length(cells)
  }

  allele_freq
}

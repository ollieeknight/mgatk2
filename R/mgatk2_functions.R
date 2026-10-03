# Helpers for mgatk2 HDF5 output.

library(hdf5r)
library(Matrix)
library(dplyr)
library(tibble)

BASES <- c("A", "C", "G", "T")

# Load one run. Matrices come back sparse, as cells x positions.
read_mgatk_hdf5 <- function(dir) {
  counts <- H5File$new(file.path(dir, "output", "counts.h5"), "r")
  meta <- H5File$new(file.path(dir, "output", "metadata.h5"), "r")
  on.exit({
    counts$close_all()
    meta$close_all()
  })

  barcodes <- counts[["barcode"]][]
  read_matrix <- function(file, name) {
    m <- as(file[[name]][, ], "CsparseMatrix")
    dimnames(m) <- list(barcodes, NULL)
    m
  }

  datasets <- c(paste0(rep(BASES, each = 2), c("_fwd", "_rev")), "tn5_cuts_fwd", "tn5_cuts_rev")
  mats <- lapply(setNames(datasets, datasets), function(n) read_matrix(counts, n))
  mats$coverage <- read_matrix(meta, "coverage")

  cells <- tibble(
    barcode = barcodes,
    mean_depth = meta[["mean_depth"]][],
    median_depth = meta[["median_depth"]][],
    max_depth = meta[["max_depth"]][],
    coverage_breadth = meta[["genome_coverage"]][],
    total_bases = meta[["total_bases"]][]
  )
  if ("barcode_metadata" %in% names(meta)) {
    group <- meta[["barcode_metadata"]]
    for (field in setdiff(names(group), "barcode")) cells[[field]] <- group[[field]][]
  }

  list(mats = mats, cells = cells, refallele = meta[["reference"]][])
}

# Combine runs that share a reference, e.g. several samples.
rbind_mgatk_data <- function(...) {
  runs <- list(...)
  mats <- lapply(setNames(names(runs[[1]]$mats), names(runs[[1]]$mats)), function(n) {
    do.call(rbind, lapply(runs, function(x) x$mats[[n]]))
  })
  list(mats = mats, cells = bind_rows(lapply(runs, `[[`, "cells")), refallele = runs[[1]]$refallele)
}

subset_cells <- function(data, keep) {
  data$mats <- lapply(data$mats, function(m) m[keep, , drop = FALSE])
  data$cells <- data$cells[keep, ]
  data
}

subset_mgatk_barcodes <- function(data, barcodes) {
  subset_cells(data, data$cells$barcode %in% barcodes)
}

filter_cells_by_coverage <- function(data, min_mean_depth, min_coverage_breadth) {
  subset_cells(
    data,
    data$cells$mean_depth >= min_mean_depth & data$cells$coverage_breadth >= min_coverage_breadth
  )
}

# The reference allele is inferred from aggregate counts, so refresh it after
# subsetting cells.
recompute_reference_alleles <- function(data) {
  totals <- sapply(BASES, function(b) {
    colSums(data$mats[[paste0(b, "_fwd")]]) + colSums(data$mats[[paste0(b, "_rev")]])
  })
  data$refallele <- ifelse(rowSums(totals) > 0, BASES[max.col(totals, "first")], "N")
  data
}

# Per-position depth, strand balance, and Tn5 summaries.
position_stats <- function(data) {
  cov <- data$mats$coverage
  fwd <- Reduce(`+`, lapply(BASES, function(b) colSums(data$mats[[paste0(b, "_fwd")]])))
  rev <- Reduce(`+`, lapply(BASES, function(b) colSums(data$mats[[paste0(b, "_rev")]])))
  tn5_fwd <- colSums(data$mats$tn5_cuts_fwd)
  tn5_rev <- colSums(data$mats$tn5_cuts_rev)
  mean_depth <- colMeans(cov)

  tibble(
    position = seq_along(mean_depth),
    mean_depth = mean_depth,
    depth_cv = sqrt(pmax(colMeans(cov^2) - mean_depth^2, 0)) / mean_depth,
    cells_covered = colSums(cov > 0),
    dropout = 1 - colSums(cov > 0) / nrow(cov),
    strand_bias = abs(fwd - rev) / pmax(fwd + rev, 1),
    tn5_fwd = tn5_fwd,
    tn5_rev = tn5_rev
  )
}

# mgatk/Signac variant statistics. VMR is the variance of per-cell allele
# frequency over the bulk frequency (plus 1e-11); with stabilise_variance, cells
# below low_coverage_threshold are held at the bulk frequency. Strand
# correlation is taken over cells carrying the allele on either strand. For
# RNA, where strand carries no information, pass min_strand_cor = -1.
identify_variants <- function(data, min_cells = 0, min_strand_cor = 0, min_vmr = 0,
                              stabilise_variance = TRUE, low_coverage_threshold = 10) {
  n_cells <- nrow(data$cells)
  cov_all <- data$mats$coverage

  per_base <- lapply(BASES, function(b) {
    idx <- which(data$refallele != b)
    fwd <- data$mats[[paste0(b, "_fwd")]][, idx, drop = FALSE]
    rev <- data$mats[[paste0(b, "_rev")]][, idx, drop = FALSE]
    alt <- fwd + rev
    cov <- cov_all[, idx, drop = FALSE]

    bulk <- colSums(alt) / colSums(cov)
    bulk[!is.finite(bulk)] <- 0

    # Allele frequency is only non-zero where alt > 0, so keep it sparse.
    af <- as(alt, "TsparseMatrix")
    af@x <- af@x / cov[cbind(af@i + 1, af@j + 1)]
    af <- as(af, "CsparseMatrix")
    if (stabilise_variance) {
      high <- cov >= low_coverage_threshold
      af <- af * high
      n_low <- n_cells - colSums(high)
      total <- colSums(af) + n_low * bulk
      total_sq <- colSums(af^2) + n_low * bulk^2
    } else {
      total <- colSums(af)
      total_sq <- colSums(af^2)
    }
    variance <- (total_sq - total^2 / n_cells) / (n_cells - 1)

    n <- colSums(alt > 0)
    f <- colSums(fwd)
    r <- colSums(rev)
    covar <- n * colSums(fwd * rev) - f * r
    spread <- (n * colSums(fwd^2) - f^2) * (n * colSums(rev^2) - r^2)
    strand_cor <- ifelse(n > 1 & spread > 0, covar / sqrt(spread), 0)

    tibble(
      position = idx,
      nucleotide = paste0(data$refallele[idx], ">", b),
      variant = paste0(idx, data$refallele[idx], ">", b),
      mean = bulk,
      vmr = variance / (bulk + 1e-11),
      n_cells_detected = n,
      n_cells_conf_detected = colSums(fwd >= 2 & rev >= 2),
      strand_correlation = strand_cor
    )
  })

  bind_rows(per_base) %>%
    filter(
      is.finite(vmr), n_cells_conf_detected >= min_cells,
      strand_correlation >= min_strand_cor, vmr > min_vmr
    ) %>%
    arrange(desc(vmr))
}

# Sparse variants x cells allele-frequency matrix.
calculate_allele_freq <- function(data, variants) {
  alt <- sub(".*>", "", variants$variant)
  counts <- vapply(seq_len(nrow(variants)), function(i) {
    p <- variants$position[i]
    data$mats[[paste0(alt[i], "_fwd")]][, p] + data$mats[[paste0(alt[i], "_rev")]][, p]
  }, numeric(nrow(data$cells)))
  af <- counts / as.matrix(data$mats$coverage[, variants$position, drop = FALSE])
  af[!is.finite(af)] <- 0
  Matrix(t(af), sparse = TRUE, dimnames = list(gsub(">", "-", variants$variant), data$cells$barcode))
}

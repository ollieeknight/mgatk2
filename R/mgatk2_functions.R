# Analysis helpers for mgatk2 HDF5 output.
# Rewritten for memory mode: arrays are fully loaded into sparse matrices.

library(hdf5r)
library(Matrix)
library(matrixStats)
library(dplyr)
library(tibble)
library(stringr)

read_mgatk_hdf5 <- function(mgatk_output_dir) {
  counts <- H5File$new(file.path(mgatk_output_dir, "output", "counts.h5"), "r")
  meta <- H5File$new(file.path(mgatk_output_dir, "output", "metadata.h5"), "r")
  
  barcodes <- counts[["barcode"]][]
  positions <- seq_len(meta[["coverage"]]$dims[2])
  
  mats <- list(coverage = as(meta[["coverage"]][,], "sparseMatrix"))
  for (b in c("A", "C", "G", "T")) {
    mats[[paste0(b, "_fwd")]] <- as(counts[[paste0(b, "_fwd")]][,], "sparseMatrix")
    mats[[paste0(b, "_rev")]] <- as(counts[[paste0(b, "_rev")]][,], "sparseMatrix")
  }
  mats$tn5_cuts_fwd <- as(counts[["tn5_cuts_fwd"]][,], "sparseMatrix")
  mats$tn5_cuts_rev <- as(counts[["tn5_cuts_rev"]][,], "sparseMatrix")
  
  for (n in names(mats)) dimnames(mats[[n]]) <- list(barcodes, positions)
  
  barcode_metadata <- NULL
  if ("barcode_metadata" %in% names(meta)) {
    barcode_metadata <- tibble(barcode = barcodes)
    for (field in setdiff(names(meta[["barcode_metadata"]]), "barcode")) {
      barcode_metadata[[field]] <- meta[["barcode_metadata"]][[field]][]
    }
  }
  
  res <- list(
    mats = mats,
    barcodes = barcodes,
    positions = positions,
    refallele = meta[["reference"]][],
    mean_depth = meta[["mean_depth"]][],
    median_depth = meta[["median_depth"]][],
    max_depth = meta[["max_depth"]][],
    genome_coverage = meta[["genome_coverage"]][],
    total_bases = meta[["total_bases"]][],
    barcode_metadata = barcode_metadata
  )
  
  counts$close_all()
  meta$close_all()
  res
}

rbind_mgatk_data <- function(...) {
  lst <- list(...)
  mats <- list()
  for (n in names(lst[[1]]$mats)) {
    mats[[n]] <- do.call(rbind, lapply(lst, function(x) x$mats[[n]]))
  }
  
  barcode_metadata <- NULL
  if (!is.null(lst[[1]]$barcode_metadata)) {
    barcode_metadata <- bind_rows(lapply(lst, function(x) x$barcode_metadata))
  }
  
  list(
    mats = mats,
    barcodes = unlist(lapply(lst, function(x) x$barcodes)),
    positions = lst[[1]]$positions,
    refallele = lst[[1]]$refallele,
    mean_depth = unlist(lapply(lst, function(x) x$mean_depth)),
    median_depth = unlist(lapply(lst, function(x) x$median_depth)),
    max_depth = unlist(lapply(lst, function(x) x$max_depth)),
    genome_coverage = unlist(lapply(lst, function(x) x$genome_coverage)),
    total_bases = unlist(lapply(lst, function(x) x$total_bases)),
    barcode_metadata = barcode_metadata
  )
}

subset_mgatk_barcodes <- function(mgatk_data, barcodes) {
  idx <- match(barcodes, mgatk_data$barcodes)
  idx <- idx[!is.na(idx)]
  
  for (n in names(mgatk_data$mats)) {
    mgatk_data$mats[[n]] <- mgatk_data$mats[[n]][idx, , drop = FALSE]
  }
  
  for (field in c("barcodes", "mean_depth", "median_depth", "max_depth", "genome_coverage", "total_bases")) {
    mgatk_data[[field]] <- mgatk_data[[field]][idx]
  }
  if (!is.null(mgatk_data$barcode_metadata)) {
    mgatk_data$barcode_metadata <- mgatk_data$barcode_metadata[idx, ]
  }
  
  mgatk_data
}

filter_cells_by_coverage <- function(mgatk_data, min_mean_coverage = 10, min_coverage_breadth = 0.5, cell_coverage_stats = NULL) {
  if (is.null(cell_coverage_stats)) {
    cell_coverage_stats <- calculate_cell_coverage_stats(mgatk_data)
  }
  
  cells_pass <- cell_coverage_stats %>%
    filter(mean_coverage >= min_mean_coverage,
           coverage_breadth >= min_coverage_breadth) %>%
    pull(barcode)
    
  subset_mgatk_barcodes(mgatk_data, cells_pass)
}

calculate_cell_coverage_stats <- function(mgatk_data, deep = FALSE) {
  stats <- tibble(
    barcode = mgatk_data$barcodes,
    mean_coverage = mgatk_data$mean_depth,
    median_coverage = mgatk_data$median_depth,
    max_coverage = mgatk_data$max_depth,
    coverage_breadth = mgatk_data$genome_coverage
  )
  if (deep) {
    cov <- as.matrix(mgatk_data$mats$coverage)
    stats$coverage_cv <- rowVars(cov) / rowMeans(cov)
    stats$positions_10x <- rowSums(cov >= 10)
    stats$positions_50x <- rowSums(cov >= 50)
  }
  stats
}

calculate_position_coverage_stats <- function(mgatk_data) {
  cov <- as.matrix(mgatk_data$mats$coverage)
  mean_cov <- colMeans(cov)
  tibble(
    position = mgatk_data$positions,
    mean_coverage = mean_cov,
    max_coverage = colMaxs(cov),
    coverage_cv = sqrt(pmax(colVars(cov), 0)) / mean_cov,
    cells_covered = colSums(cov > 0),
    cells_10x = colSums(cov >= 10),
    cells_50x = colSums(cov >= 50),
    coverage_dropout = 1 - colSums(cov > 0) / nrow(cov)
  )
}

recompute_reference_alleles <- function(mgatk_data) {
  bases <- c("A", "C", "G", "T")
  old_ref <- mgatk_data$refallele
  
  # Vectorized calculation of base totals
  base_counts <- lapply(bases, function(b) {
    Matrix::colSums(mgatk_data$mats[[paste0(b, "_fwd")]] + mgatk_data$mats[[paste0(b, "_rev")]])
  })
  # Bind into a matrix (positions x 4)
  base_mat <- do.call(cbind, base_counts)
  
  # Find max base for each position
  max_idx <- max.col(base_mat, ties.method = "first")
  new_ref <- bases[max_idx]
  
  n_changed <- sum(old_ref != new_ref)
  message("Updated ", n_changed, " reference alleles")
  
  mgatk_data$refallele <- new_ref
  mgatk_data
}

calculate_strand_coverage_stats <- function(mgatk_data) {
  n_cells <- length(mgatk_data$barcodes)
  fwd <- numeric(length(mgatk_data$positions))
  rev <- numeric(length(mgatk_data$positions))
  for (b in c("A", "C", "G", "T")) {
    fwd <- fwd + Matrix::colSums(mgatk_data$mats[[paste0(b, "_fwd")]])
    rev <- rev + Matrix::colSums(mgatk_data$mats[[paste0(b, "_rev")]])
  }
  total <- Matrix::colSums(mgatk_data$mats$coverage)
  
  tibble(
    position = mgatk_data$positions,
    mean_fwd_coverage = fwd / n_cells,
    mean_rev_coverage = rev / n_cells,
    total_coverage = total / n_cells,
    strand_balance = (fwd / n_cells) / pmax(total / n_cells, 1),
    strand_imbalance = abs((fwd - rev) / n_cells)
  )
}

calculate_transposition_stats <- function(mgatk_data) {
  n_cells <- length(mgatk_data$barcodes)
  fwd <- Matrix::colSums(mgatk_data$mats$tn5_cuts_fwd)
  rev <- Matrix::colSums(mgatk_data$mats$tn5_cuts_rev)
  total <- Matrix::colSums(mgatk_data$mats$coverage)
  
  tibble(
    position = mgatk_data$positions,
    tn5_cuts_total = fwd + rev,
    tn5_cuts_fwd = fwd,
    tn5_cuts_rev = rev,
    total_bases = total,
    tn5_frequency = ifelse(total > 0, (fwd + rev) / total, 0),
    cells_with_tn5_cuts = Matrix::colSums(mgatk_data$mats$tn5_cuts_fwd + mgatk_data$mats$tn5_cuts_rev > 0),
    mean_tn5_cuts_per_cell = (fwd + rev) / n_cells,
    strand_bias = abs(fwd - rev) / pmax(fwd + rev, 1)
  )
}

identify_variants <- function(mgatk_data, min_cells = 0, min_strand_cor = 0, min_vmr = 0,
                              stabilise_variance = FALSE, low_coverage_threshold = 10) {
  bases <- c("A", "C", "G", "T")
  n_cells <- length(mgatk_data$barcodes)
  cov <- mgatk_data$mats$coverage
  cov_total <- Matrix::colSums(cov)
  
  results <- lapply(bases, function(b) {
    # Match the old behavior perfectly: only exclude `b`, don't exclude "N" explicitly
    idx <- which(mgatk_data$refallele != b)
    if (length(idx) == 0) return(NULL)
    
    fwd <- mgatk_data$mats[[paste0(b, "_fwd")]][, idx, drop = FALSE]
    rev <- mgatk_data$mats[[paste0(b, "_rev")]][, idx, drop = FALSE]
    alt <- fwd + rev
    
    bulk <- Matrix::colSums(alt) / cov_total[idx]
    bulk[is.na(bulk)] <- 0
    
    if (stabilise_variance) {
      low <- cov[, idx, drop = FALSE] < low_coverage_threshold
      af <- as.matrix(alt / cov[, idx, drop = FALSE])
      af[is.na(af)] <- 0
      bulk_mat <- matrix(rep(bulk, n_cells), nrow = n_cells, byrow = TRUE)
      af[as.matrix(low)] <- bulk_mat[as.matrix(low)]
      variance <- colVars(af)
    } else {
      af <- as.matrix(alt / cov[, idx, drop = FALSE])
      af[is.na(af)] <- 0
      variance <- colVars(af)
    }
    
    # Old calculation divided variance / bulk, meaning no 1e-11 pseudocount.
    # Reverting to exact old math to guarantee matching VMRs
    vmr <- variance / bulk
    
    # Match old strand correlation perfectly: cor(fwd, rev) over cells where fwd > 0
    seen <- (fwd > 0)
    n <- Matrix::colSums(seen)
    f <- Matrix::colSums(fwd * seen)
    r <- Matrix::colSums(rev * seen)
    
    covar <- n * Matrix::colSums(fwd * rev * seen) - f * r
    spread <- (n * Matrix::colSums(fwd^2 * seen) - f^2) * (n * Matrix::colSums(rev^2 * seen) - r^2)
    strand_cor <- ifelse(n > 1 & spread > 0, covar / sqrt(spread), 0)
    strand_cor[is.na(strand_cor)] <- 0
    
    tibble(
      position = mgatk_data$positions[idx],
      nucleotide = paste0(mgatk_data$refallele[idx], ">", b),
      variant = paste0(mgatk_data$positions[idx], mgatk_data$refallele[idx], ">", b),
      mean = bulk,
      vmr = vmr,
      n_cells_detected = Matrix::colSums(alt > 0),
      n_cells_conf_detected = Matrix::colSums(fwd >= 2 & rev >= 2),
      strand_correlation = strand_cor
    )
  })
  
  bind_rows(results) %>%
    filter(is.finite(vmr), n_cells_conf_detected >= min_cells, strand_correlation >= min_strand_cor, vmr > min_vmr) %>%
    arrange(desc(vmr))
}

calculate_allele_freq <- function(mgatk_data, variants, min_coverage = 1) {
  if (nrow(variants) == 0) {
    mat <- Matrix::Matrix(0, nrow = 0, ncol = length(mgatk_data$barcodes), sparse = TRUE)
    colnames(mat) <- mgatk_data$barcodes
    return(mat)
  }
  
  variant_names <- gsub(">", "-", variants$variant)
  
  allele_freq <- Matrix::Matrix(0, nrow = nrow(variants), ncol = length(mgatk_data$barcodes), 
                                sparse = TRUE)
  rownames(allele_freq) <- variant_names
  colnames(allele_freq) <- mgatk_data$barcodes
  
  for (i in seq_len(nrow(variants))) {
    parts <- str_match(variants$variant[i], "(\\d+)([ACGT])>([ACGT])")
    pos <- as.numeric(parts[2])
    alt_base <- parts[4]
    
    pos_idx <- which(mgatk_data$positions == pos)
    if (length(pos_idx) == 0) next
    
    fwd <- mgatk_data$mats[[paste0(alt_base, "_fwd")]][, pos_idx]
    rev <- mgatk_data$mats[[paste0(alt_base, "_rev")]][, pos_idx]
    cov <- mgatk_data$mats$coverage[, pos_idx]
    
    af <- ifelse(cov >= min_coverage, (fwd + rev) / cov, 0)
    allele_freq[i, ] <- af
  }
  
  allele_freq
}

probe_variant <- function(mgatk_data, variant_str, stabilise_variance = FALSE, low_coverage_threshold = 10) {
  # Parse variant_str, e.g., "13710A>G"
  matches <- regmatches(variant_str, regexec("^(\\d+)([ACGTN])>([ACGTN])$", variant_str))[[1]]
  if (length(matches) != 4) {
    stop("variant_str must be in format 'POSRef>Alt', e.g. '13710A>G'")
  }
  
  pos <- as.integer(matches[2])
  ref <- matches[3]
  alt_base <- matches[4]
  
  idx <- which(mgatk_data$positions == pos)
  if (length(idx) == 0) {
    stop("Position ", pos, " not found in mgatk_data")
  }
  
  if (mgatk_data$refallele[idx] != ref) {
    warning("Reference allele at position ", pos, " is ", mgatk_data$refallele[idx], ", but variant_str specifies ", ref)
  }
  
  cov <- as.vector(mgatk_data$mats$coverage[, idx])
  cov_total <- sum(cov)
  
  fwd <- as.vector(mgatk_data$mats[[paste0(alt_base, "_fwd")]][, idx])
  rev <- as.vector(mgatk_data$mats[[paste0(alt_base, "_rev")]][, idx])
  alt <- fwd + rev
  
  bulk <- sum(alt) / cov_total
  if (is.na(bulk)) bulk <- 0
  
  if (stabilise_variance) {
    low <- cov < low_coverage_threshold
    af <- alt / cov
    af[is.na(af)] <- 0
    af[low] <- bulk
    variance <- var(af)
  } else {
    af <- alt / cov
    af[is.na(af)] <- 0
    variance <- var(af)
  }
  
  vmr <- variance / bulk
  
  seen <- (fwd > 0)
  n <- sum(seen)
  f <- sum(fwd[seen])
  r <- sum(rev[seen])
  
  covar <- n * sum(fwd[seen] * rev[seen]) - f * r
  spread <- (n * sum(fwd[seen]^2) - f^2) * (n * sum(rev[seen]^2) - r^2)
  strand_cor <- ifelse(n > 1 & spread > 0, covar / sqrt(spread), 0)
  
  stats <- tibble::tibble(
    position = pos,
    nucleotide = paste0(ref, ">", alt_base),
    variant = variant_str,
    mean = bulk,
    vmr = vmr,
    n_cells_detected = sum(alt > 0),
    n_cells_conf_detected = sum(fwd >= 2 & rev >= 2),
    strand_correlation = strand_cor
  )
  
  list(
    stats = stats,
    cell_allele_freq = setNames(af, mgatk_data$barcodes),
    cell_coverage = setNames(cov, mgatk_data$barcodes),
    cell_alt_counts = setNames(alt, mgatk_data$barcodes),
    cell_fwd_counts = setNames(fwd, mgatk_data$barcodes),
    cell_rev_counts = setNames(rev, mgatk_data$barcodes)
  )
}

identify_variants_v2 <- function(mgatk_data, mode = c("ATAC", "RNA"), min_cells = 0, min_strand_cor = NULL, min_vmr = 0.01,
                                 stabilise_variance = TRUE, low_coverage_threshold = 10) {
  mode <- match.arg(mode)
  
  if (is.null(min_strand_cor)) {
    min_strand_cor <- ifelse(mode == "ATAC", 0.65, -1)
  }
  
  bases <- c("A", "C", "G", "T")
  n_cells <- length(mgatk_data$barcodes)
  cov <- mgatk_data$mats$coverage
  cov_total <- Matrix::colSums(cov)
  
  results <- lapply(bases, function(b) {
    idx <- which(mgatk_data$refallele != b)
    if (length(idx) == 0) return(NULL)
    
    # Subset to the current base
    fwd <- mgatk_data$mats[[paste0(b, "_fwd")]][, idx, drop = FALSE]
    rev <- mgatk_data$mats[[paste0(b, "_rev")]][, idx, drop = FALSE]
    alt <- fwd + rev
    cov_sub <- cov[, idx, drop = FALSE]
    
    # Bulk AF
    bulk <- Matrix::colSums(alt) / cov_total[idx]
    bulk[is.na(bulk)] <- 0
    
    # Sparse allele frequency (only computes where alt > 0)
    # This prevents the memory bomb of as.matrix()
    af <- alt
    alt_trip <- summary(alt)
    if (nrow(alt_trip) > 0) {
      cov_vals <- cov_sub[cbind(alt_trip$i, alt_trip$j)]
      af_x <- alt_trip$x / cov_vals
      af <- Matrix::sparseMatrix(i = alt_trip$i, j = alt_trip$j, x = af_x, dims = dim(alt))
    }
    
    if (stabilise_variance) {
      # High coverage cells (sparse logical matrix)
      high_cov <- cov_sub >= low_coverage_threshold
      N_H <- Matrix::colSums(high_cov)
      N_L <- n_cells - N_H
      
      af_high <- af * high_cov
      sum_X_H <- Matrix::colSums(af_high)
      
      af_high_sq <- af_high
      af_high_sq@x <- af_high_sq@x^2
      sum_X2_H <- Matrix::colSums(af_high_sq)
      
      S <- sum_X_H + N_L * bulk
      S2 <- sum_X2_H + N_L * (bulk^2)
      
      variance <- (S2 - (S^2 / n_cells)) / (n_cells - 1)
    } else {
      sum_X <- Matrix::colSums(af)
      af_sq <- af
      af_sq@x <- af_sq@x^2
      sum_X2 <- Matrix::colSums(af_sq)
      variance <- (sum_X2 - (sum_X^2 / n_cells)) / (n_cells - 1)
    }
    
    # Pseudocount protects against div-by-zero
    vmr <- variance / (bulk + 1e-11)
    
    # Strand correlation (only relevant if not RNA)
    if (mode == "ATAC") {
      # Compute correlation across cells where variant is detected (alt > 0)
      # This fixes the biological flaw of dropping cells with no fwd reads
      seen <- (alt > 0)
      n <- Matrix::colSums(seen)
      
      f <- Matrix::colSums(fwd)
      r <- Matrix::colSums(rev)
      
      f2 <- fwd
      f2@x <- f2@x^2
      r2 <- rev
      r2@x <- r2@x^2
      
      sum_f2 <- Matrix::colSums(f2)
      sum_r2 <- Matrix::colSums(r2)
      
      covar <- n * Matrix::colSums(fwd * rev) - f * r
      spread <- (n * sum_f2 - f^2) * (n * sum_r2 - r^2)
      
      strand_cor <- ifelse(n > 1 & spread > 0, covar / sqrt(spread), 0)
      strand_cor[is.na(strand_cor)] <- 0
    } else {
      strand_cor <- NA_real_
    }
    
    tibble::tibble(
      position = mgatk_data$positions[idx],
      nucleotide = paste0(mgatk_data$refallele[idx], ">", b),
      variant = paste0(mgatk_data$positions[idx], mgatk_data$refallele[idx], ">", b),
      mean = bulk,
      vmr = vmr,
      n_cells_detected = Matrix::colSums(alt > 0),
      n_cells_conf_detected = Matrix::colSums(fwd >= 2 & rev >= 2),
      strand_correlation = strand_cor
    )
  })
  
  res <- dplyr::bind_rows(results) %>%
    dplyr::filter(is.finite(vmr), n_cells_conf_detected >= min_cells, vmr > min_vmr)
    
  if (mode == "ATAC") {
    res <- res %>% dplyr::filter(strand_correlation >= min_strand_cor)
  }
  
  res %>% dplyr::arrange(dplyr::desc(vmr))
}

read_mgatk_hdf5 <- function(mgatk_output_dir) {
  
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required. Install with: install.packages('hdf5r')")
  }

  counts_file <- file.path(mgatk_output_dir, "output", "counts.h5")
  metadata_file <- file.path(mgatk_output_dir, "output", "metadata.h5")
  
  if (!file.exists(counts_file)) {
    stop("Counts file not found: ", counts_file)
  }
  
  if (!file.exists(metadata_file)) {
    stop("Metadata file not found: ", metadata_file)
  }
  
  counts_h5 <- hdf5r::H5File$new(counts_file, mode = "r")
  metadata_h5 <- hdf5r::H5File$new(metadata_file, mode = "r")
  
  barcodes <- counts_h5[["barcode"]][]
  
  message("Reading count matrices...")
  A_fwd_data <- counts_h5[["A_fwd"]][,]
  n_cells <- nrow(A_fwd_data)
  n_positions <- ncol(A_fwd_data)
  
  message("  Found ", n_cells, " cells across ", n_positions, " positions")
  
  positions <- 1:n_positions
  
  counts <- list(
    A_fwd = A_fwd_data,
    A_rev = counts_h5[["A_rev"]][,],
    C_fwd = counts_h5[["C_fwd"]][,],
    C_rev = counts_h5[["C_rev"]][,],
    G_fwd = counts_h5[["G_fwd"]][,],
    G_rev = counts_h5[["G_rev"]][,],
    T_fwd = counts_h5[["T_fwd"]][,],
    T_rev = counts_h5[["T_rev"]][,],
    tn5_cuts_fwd = counts_h5[["tn5_cuts_fwd"]][,],
    tn5_cuts_rev = counts_h5[["tn5_cuts_rev"]][,],
    coverage = metadata_h5[["coverage"]][,]
  )
  
  mean_depth <- metadata_h5[["mean_depth"]][]
  names(mean_depth) <- barcodes
  
  median_depth <- metadata_h5[["median_depth"]][]
  names(median_depth) <- barcodes
  
  max_depth <- metadata_h5[["max_depth"]][]
  names(max_depth) <- barcodes
  
  genome_coverage <- metadata_h5[["genome_coverage"]][]
  names(genome_coverage) <- barcodes
  
  total_bases <- metadata_h5[["total_bases"]][]
  names(total_bases) <- barcodes
  
  refallele <- metadata_h5[["reference"]][]
  
  barcode_metadata <- NULL
  if ("barcode_metadata" %in% names(metadata_h5)) {
    bm <- metadata_h5[["barcode_metadata"]]
    
    metadata_fields <- names(bm)
    barcode_metadata <- tibble(barcode = barcodes)
    
    for (field in metadata_fields) {
      if (field != "barcode") {
        tryCatch({
          barcode_metadata[[field]] <- bm[[field]][]
        }, error = function(e) {
          warning(paste("Could not read barcode_metadata field:", field))
        })
      }
    }
  }
  
  counts_h5$close_all()
  metadata_h5$close_all()
  
  list(
    counts = counts,
    mean_depth = mean_depth,
    median_depth = median_depth,
    max_depth = max_depth,
    genome_coverage = genome_coverage,
    total_bases = total_bases,
    refallele = refallele,
    positions = positions,
    barcodes = barcodes,
    barcode_metadata = barcode_metadata
  )
}

subset_mgatk_barcodes <- function(mgatk_data, barcodes) {
  
  keep_idx <- which(mgatk_data$barcodes %in% barcodes)
  
  if (length(keep_idx) == 0) {
    stop("No matching barcodes found")
  }
  
  message("  Found ", length(keep_idx), " matching barcodes (", 
          round(100 * length(keep_idx) / length(mgatk_data$barcodes), 1), "%)")
  
  counts_subset <- lapply(mgatk_data$counts, function(mat) {
    if (!is.null(mat)) {
      mat[keep_idx, , drop = FALSE]
    } else {
      NULL
    }
  })
  
  barcode_metadata_subset <- NULL
  if (!is.null(mgatk_data$barcode_metadata)) {
    barcode_metadata_subset <- mgatk_data$barcode_metadata[keep_idx, ]
  }
  
  result <- list(
    counts = counts_subset,
    mean_depth = mgatk_data$mean_depth[keep_idx],
    median_depth = mgatk_data$median_depth[keep_idx],
    max_depth = mgatk_data$max_depth[keep_idx],
    genome_coverage = mgatk_data$genome_coverage[keep_idx],
    total_bases = mgatk_data$total_bases[keep_idx],
    refallele = mgatk_data$refallele,
    positions = mgatk_data$positions,
    barcodes = mgatk_data$barcodes[keep_idx],
    barcode_metadata = barcode_metadata_subset
  )
  
  return(result)
}

filter_cells_by_coverage <- function(mgatk_data, 
                                      min_mean_coverage = 10, 
                                      min_coverage_breadth = 0.5,
                                      cell_coverage_stats = NULL) {
  
  if (is.null(cell_coverage_stats)) {
    cell_coverage_stats <- calculate_cell_coverage_stats(mgatk_data)
  }
  
  cells_pass <- cell_coverage_stats %>%
    filter(mean_coverage >= min_mean_coverage,
           coverage_breadth >= min_coverage_breadth) %>%
    pull(barcode)
  
  keep_idx <- which(mgatk_data$barcodes %in% cells_pass)
  
  counts_filtered <- lapply(mgatk_data$counts, function(mat) {
    if (!is.null(mat)) mat[keep_idx, , drop = FALSE] else NULL
  })
  
  barcode_metadata_filtered <- NULL
  if (!is.null(mgatk_data$barcode_metadata)) {
    barcode_metadata_filtered <- mgatk_data$barcode_metadata[keep_idx, ]
  }
  
  message("Kept ", length(keep_idx), " / ", length(mgatk_data$barcodes), " cells")
  
  list(
    counts = counts_filtered,
    mean_depth = mgatk_data$mean_depth[keep_idx],
    median_depth = mgatk_data$median_depth[keep_idx],
    max_depth = mgatk_data$max_depth[keep_idx],
    genome_coverage = mgatk_data$genome_coverage[keep_idx],
    total_bases = mgatk_data$total_bases[keep_idx],
    refallele = mgatk_data$refallele,
    positions = mgatk_data$positions,
    barcodes = mgatk_data$barcodes[keep_idx],
    barcode_metadata = barcode_metadata_filtered
  )
}

calculate_cell_coverage_stats <- function(mgatk_data) {
  coverage_matrix <- mgatk_data$counts$coverage
  
  cell_stats <- tibble(
    barcode = mgatk_data$barcodes,
    mean_coverage = rowMeans(coverage_matrix),
    median_coverage = apply(coverage_matrix, 1, median),
    max_coverage = apply(coverage_matrix, 1, max),
    coverage_cv = apply(coverage_matrix, 1, function(x) sd(x) / mean(x)),
    positions_covered = apply(coverage_matrix > 0, 1, sum),
    positions_10x = apply(coverage_matrix >= 10, 1, sum),
    positions_50x = apply(coverage_matrix >= 50, 1, sum),
    coverage_breadth = apply(coverage_matrix > 0, 1, sum) / ncol(coverage_matrix)
  )
  
  return(cell_stats)
}

calculate_position_coverage_stats <- function(mgatk_data) {
  coverage_matrix <- mgatk_data$counts$coverage
  
  position_stats <- tibble(
    position = mgatk_data$positions,
    mean_coverage = colMeans(coverage_matrix),
    median_coverage = apply(coverage_matrix, 2, median),
    max_coverage = apply(coverage_matrix, 2, max),
    coverage_cv = apply(coverage_matrix, 2, function(x) sd(x) / mean(x)),
    cells_covered = apply(coverage_matrix > 0, 2, sum),
    cells_10x = apply(coverage_matrix >= 10, 2, sum),
    cells_50x = apply(coverage_matrix >= 50, 2, sum),
    coverage_dropout = apply(coverage_matrix == 0, 2, sum) / nrow(coverage_matrix)
  )
  
  return(position_stats)
}

recompute_reference_alleles <- function(mgatk_data) {
  
  bases <- c("A", "C", "G", "T")
  old_ref <- mgatk_data$refallele
  new_ref <- character(length(mgatk_data$positions))
  
  for (i in seq_along(mgatk_data$positions)) {
    base_totals <- sapply(bases, function(base) {
      sum(mgatk_data$counts[[paste0(base, "_fwd")]][, i] + 
          mgatk_data$counts[[paste0(base, "_rev")]][, i])
    })
    new_ref[i] <- bases[which.max(base_totals)]
  }
  
  n_changed <- sum(old_ref != new_ref)
  message("Updated ", n_changed, " reference alleles")
  
  mgatk_data$refallele <- new_ref
  mgatk_data
}

calculate_strand_coverage_stats <- function(mgatk_data) {
  
  fwd_coverage <- mgatk_data$counts$A_fwd + mgatk_data$counts$C_fwd + 
                  mgatk_data$counts$G_fwd + mgatk_data$counts$T_fwd
  
  rev_coverage <- mgatk_data$counts$A_rev + mgatk_data$counts$C_rev + 
                  mgatk_data$counts$G_rev + mgatk_data$counts$T_rev
  
  tibble(
    position = mgatk_data$positions,
    mean_fwd_coverage = colMeans(fwd_coverage),
    mean_rev_coverage = colMeans(rev_coverage),
    total_coverage = colMeans(mgatk_data$counts$coverage)
  ) %>%
    mutate(
      strand_balance = mean_fwd_coverage / pmax(total_coverage, 1),
      strand_imbalance = abs(mean_fwd_coverage - mean_rev_coverage)
    )
}

calculate_transposition_stats <- function(mgatk_data) {
  
  if (is.null(mgatk_data$counts$tn5_cuts_fwd) || is.null(mgatk_data$counts$tn5_cuts_rev)) {
    return(tibble())
  }
  
  tn5_fwd <- mgatk_data$counts$tn5_cuts_fwd
  tn5_rev <- mgatk_data$counts$tn5_cuts_rev
  tn5_total <- tn5_fwd + tn5_rev
  cov <- mgatk_data$counts$coverage
  
  total_bases <- colSums(cov)
  
  tibble(
    position = mgatk_data$positions,
    tn5_cuts_total = colSums(tn5_total),
    tn5_cuts_fwd = colSums(tn5_fwd),
    tn5_cuts_rev = colSums(tn5_rev),
    total_bases = total_bases,
    tn5_frequency = ifelse(total_bases > 0, colSums(tn5_total) / total_bases, 0),
    cells_with_tn5_cuts = colSums(tn5_total > 0),
    mean_tn5_cuts_per_cell = colMeans(tn5_total),
    strand_bias = abs(tn5_cuts_fwd - tn5_cuts_rev) / pmax(tn5_cuts_fwd + tn5_cuts_rev, 1)
  )
}

identify_variants <- function(mgatk_data, min_cells = 0, min_strand_cor = 0, min_vmr = 0, 
                              stabilize_variance = FALSE, low_coverage_threshold = 10) {
  
  bases <- c("A", "C", "G", "T")
  results <- list()
  coverage <- mgatk_data$counts$coverage
  
  for (base in bases) {
    fwd <- mgatk_data$counts[[paste0(base, "_fwd")]]
    rev <- mgatk_data$counts[[paste0(base, "_rev")]]
    
    alt_pos <- which(mgatk_data$refallele != base)
    
    for (i in alt_pos) {
      fwd_counts <- fwd[, i]
      rev_counts <- rev[, i]
      cov <- coverage[, i]
      total <- fwd_counts + rev_counts
      
      cells_conf <- sum(fwd_counts >= 2 & rev_counts >= 2)
      if (cells_conf < min_cells) next
      
      bulk_af <- sum(total) / sum(cov)
      if (is.na(bulk_af)) next
      
      fwd_idx <- fwd_counts > 0
      strand_cor <- if (sum(fwd_idx) > 1) {
        cor(fwd_counts[fwd_idx], rev_counts[fwd_idx], use = "complete.obs")
      } else 0
      strand_cor[is.na(strand_cor)] <- 0
      
      af <- ifelse(cov > 0, total / cov, 0)
      if (stabilize_variance) {
        af[cov < low_coverage_threshold] <- bulk_af
      }
      vmr <- var(af) / bulk_af
      
      results[[length(results) + 1]] <- tibble(
        position = mgatk_data$positions[i],
        nucleotide = paste0(mgatk_data$refallele[i], ">", base),
        variant = paste0(mgatk_data$positions[i], mgatk_data$refallele[i], ">", base),
        mean = bulk_af,
        vmr = vmr,
        n_cells_detected = sum(total > 0),
        n_cells_conf_detected = cells_conf,
        strand_correlation = strand_cor
      )
    }
  }
  
  if (length(results) == 0) return(tibble())
  
  bind_rows(results) %>%
    filter(n_cells_conf_detected >= min_cells,
           strand_correlation >= min_strand_cor,
           vmr > min_vmr) %>%
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
    
    fwd <- mgatk_data$counts[[paste0(alt_base, "_fwd")]][, pos_idx]
    rev <- mgatk_data$counts[[paste0(alt_base, "_rev")]][, pos_idx]
    cov <- mgatk_data$counts$coverage[, pos_idx]
    
    allele_freq[i, ] <- ifelse(cov >= min_coverage, (fwd + rev) / cov, 0)
  }
  
  allele_freq
}
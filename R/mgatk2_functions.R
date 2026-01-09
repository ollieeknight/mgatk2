read_mgatk_hdf5 <- function(output_dir) {

  counts_file <- file.path(output_dir, "output", "counts.h5")
  metadata_file <- file.path(output_dir, "output", "metadata.h5")
  
  counts_h5 <- H5File$new(counts_file, mode = "r")
  metadata_h5 <- H5File$new(metadata_file, mode = "r")
  
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

filter_cells_by_coverage <- function(mgatk_data, 
                                     min_mean_coverage = 10, 
                                     min_coverage_breadth = 0.5,
                                     cell_coverage_stats = NULL) {
    
    # Calculate coverage stats if not provided (avoids redundant calculation)
    if (is.null(cell_coverage_stats)) {
        cell_coverage_stats <- calculate_cell_coverage_stats(mgatk_data)
    }
    
    # Identify cells passing filters
    cells_pass <- cell_coverage_stats %>%
        filter(mean_coverage >= min_mean_coverage,
               coverage_breadth >= min_coverage_breadth) %>%
        pull(barcode)
    
    # Get indices of passing cells
    cell_indices <- which(mgatk_data$barcodes %in% cells_pass)
    
    # Subset ALL count matrices (preserving all strand-specific and tn5 data)
    counts_filtered <- lapply(mgatk_data$counts, function(mat) {
        if (!is.null(mat)) {
            mat[cell_indices, , drop = FALSE]
        } else {
            NULL
        }
    })
    
    # Subset barcode metadata if it exists
    barcode_metadata_filtered <- NULL
    if (!is.null(mgatk_data$barcode_metadata)) {
        barcode_metadata_filtered <- mgatk_data$barcode_metadata[cell_indices, ]
    }
    
    # Preserve ALL fields from original mgatk_data structure
    mgatk_filtered <- list(
        counts = counts_filtered,
        mean_depth = mgatk_data$mean_depth[cell_indices],
        median_depth = mgatk_data$median_depth[cell_indices],
        max_depth = mgatk_data$max_depth[cell_indices],
        genome_coverage = mgatk_data$genome_coverage[cell_indices],
        total_bases = mgatk_data$total_bases[cell_indices],
        refallele = mgatk_data$refallele,
        positions = mgatk_data$positions,
        barcodes = mgatk_data$barcodes[cell_indices],
        barcode_metadata = barcode_metadata_filtered
    )
    
    message(sprintf("Filtered %d cells",
                    length(cell_indices), 
                    100 * length(cell_indices) / length(mgatk_data$barcodes),
                    min_mean_coverage, 
                    min_coverage_breadth))
    
    return(mgatk_filtered)
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

calculate_transposition_stats <- function(mgatk_data) {
  tn5_cuts_fwd <- mgatk_data$counts$tn5_cuts_fwd
  tn5_cuts_rev <- mgatk_data$counts$tn5_cuts_rev
  coverage <- mgatk_data$counts$coverage
  
  if (is.null(tn5_cuts_fwd) || is.null(tn5_cuts_rev)) {
    return(tibble(
      position = integer(),
      tn5_cuts_total = integer(),
      tn5_cuts_fwd = integer(),
      tn5_cuts_rev = integer(),
      total_bases = integer(),
      tn5_frequency = numeric(),
      cells_with_tn5_cuts = integer(),
      mean_tn5_cuts_per_cell = numeric(),
      strand_bias = numeric()
    ))
  }
  
  total_tn5_cuts <- tn5_cuts_fwd + tn5_cuts_rev
  
  transposition_stats <- tibble(
    position = mgatk_data$positions,
    tn5_cuts_total = colSums(total_tn5_cuts, na.rm = TRUE),
    tn5_cuts_fwd = colSums(tn5_cuts_fwd, na.rm = TRUE),
    tn5_cuts_rev = colSums(tn5_cuts_rev, na.rm = TRUE),
    total_bases = colSums(coverage, na.rm = TRUE),
    tn5_frequency = ifelse(colSums(coverage, na.rm = TRUE) > 0, 
                           colSums(total_tn5_cuts, na.rm = TRUE) / colSums(coverage, na.rm = TRUE), 0),
    cells_with_tn5_cuts = apply(total_tn5_cuts > 0, 2, sum),
    mean_tn5_cuts_per_cell = colMeans(total_tn5_cuts, na.rm = TRUE)
  ) %>%
    mutate(
      strand_bias = ifelse(tn5_cuts_fwd + tn5_cuts_rev > 0, 
                           abs(tn5_cuts_fwd - tn5_cuts_rev) / (tn5_cuts_fwd + tn5_cuts_rev), 0)
    )
  
  return(transposition_stats)
}

identify_variants <- function(mgatk_data, min_cells = 0, min_strand_cor = 0, min_vmr = 0, 
                              stabilize_variance = FALSE, low_coverage_threshold = 10) {
  positions <- mgatk_data$positions
  refallele <- mgatk_data$refallele
  barcodes <- mgatk_data$barcodes
  n_cells <- length(barcodes)
  n_positions <- length(positions)
  
  if (min_cells > 0) message("  Minimum cells: ", min_cells)
  if (min_strand_cor > 0) message("  Minimum strand correlation: ", min_strand_cor)
  if (min_vmr > 0) message("  Minimum VMR: ", min_vmr)
  
  bases <- c("A", "C", "G", "T")
  results_list <- list()
  
  coverage <- mgatk_data$counts$coverage
  
  for (base in bases) {
    fwd <- mgatk_data$counts[[paste0(base, "_fwd")]]
    rev <- mgatk_data$counts[[paste0(base, "_rev")]]
    
    is_ref <- refallele == base
    alt_positions <- which(!is_ref)
    
    for (pos_idx in alt_positions) {
      pos <- positions[pos_idx]
      ref_base <- refallele[pos_idx]
      
      fwd_counts <- fwd[, pos_idx]
      rev_counts <- rev[, pos_idx]
      cov <- coverage[, pos_idx]
      
      total_counts <- fwd_counts + rev_counts
      af <- ifelse(cov > 0, total_counts / cov, 0)
      
      cells_detected <- sum(total_counts > 0)
      cells_confident <- sum(fwd_counts >= 2 & rev_counts >= 2)
      
      if (cells_confident < min_cells) next
      
      bulk_af <- sum(total_counts) / sum(cov)
      if (is.na(bulk_af) || is.nan(bulk_af)) next
      
      fwd_idx <- which(fwd_counts > 0)
      if (length(fwd_idx) > 1) {
        strand_cor <- cor(fwd_counts[fwd_idx], rev_counts[fwd_idx], 
                          method = "pearson", use = "complete.obs")
      } else {
        strand_cor <- 0
      }
      if (is.na(strand_cor)) strand_cor <- 0
      
      if (stabilize_variance) {
        af_stabilized <- af
        low_cov_idx <- which(cov < low_coverage_threshold)
        af_stabilized[low_cov_idx] <- bulk_af
        var_af <- var(af_stabilized)
      } else {
        var_af <- var(af)
      }
      vmr <- var_af / bulk_af
      
      results_list[[length(results_list) + 1]] <- tibble(
        position = pos,
        nucleotide = paste0(ref_base, ">", base),
        variant = paste0(pos, ref_base, ">", base),
        mean = bulk_af,
        vmr = vmr,
        n_cells_detected = cells_detected,
        n_cells_conf_detected = cells_confident,
        strand_correlation = strand_cor
      )
    }
  }
  
  if (length(results_list) == 0) {
    return(tibble(
      position = integer(),
      nucleotide = character(),
      variant = character(),
      mean = numeric(),
      vmr = numeric(),
      n_cells_detected = integer(),
      n_cells_conf_detected = integer(),
      strand_correlation = numeric()
    ))
  }
  
  all_variants <- bind_rows(results_list)
  
  filtered_variants <- all_variants %>%
    filter(
      n_cells_conf_detected >= min_cells,
      strand_correlation >= min_strand_cor,
      vmr > min_vmr
    ) %>%
    arrange(desc(vmr))
  
  return(filtered_variants)
}

calculate_allele_freq <- function(mgatk_data, variants, min_coverage = 1) {
  
  positions <- mgatk_data$positions
  barcodes <- mgatk_data$barcodes
  n_cells <- length(barcodes)
  
  if (nrow(variants) == 0) {
    empty_matrix <- matrix(0, nrow = 0, ncol = n_cells)
    empty_matrix <- as(empty_matrix, "dgCMatrix")
    colnames(empty_matrix) <- barcodes
    return(empty_matrix)
  }
  
  
  valid_variants_idx <- c()
  for (i in seq_len(nrow(variants))) {
    var <- variants$variant[i]
    parts <- str_match(var, "(\\d+)([ACGT])>([ACGT])")
    pos <- as.numeric(parts[2])
    
    pos_idx <- which(positions == pos)
    if (length(pos_idx) > 0) {
      valid_variants_idx <- c(valid_variants_idx, i)
    }
  }
  
  if (length(valid_variants_idx) == 0) {
    empty_matrix <- matrix(0, nrow = 0, ncol = n_cells)
    empty_matrix <- as(empty_matrix, "dgCMatrix")
    colnames(empty_matrix) <- barcodes
    return(empty_matrix)
  }
  
  variants_valid <- variants[valid_variants_idx, ]
  n_valid_variants <- nrow(variants_valid)
  
  variant_names <- gsub(">", "-", variants_valid$variant)
  
  allele_freq <- matrix(0, nrow = n_valid_variants, ncol = n_cells)
  allele_freq <- as(allele_freq, "dgCMatrix")
  rownames(allele_freq) <- variant_names
  colnames(allele_freq) <- barcodes
  
  for (i in seq_len(n_valid_variants)) {
    var <- variants_valid$variant[i]
    parts <- str_match(var, "(\\d+)([ACGT])>([ACGT])")
    pos <- as.numeric(parts[2])
    alt_base <- parts[4]
    
    pos_idx <- which(positions == pos)
    
    fwd <- mgatk_data$counts[[paste0(alt_base, "_fwd")]][, pos_idx]
    rev <- mgatk_data$counts[[paste0(alt_base, "_rev")]][, pos_idx]
    cov <- mgatk_data$counts$coverage[, pos_idx]
    
    total_alt <- fwd + rev
    
    af <- ifelse(cov >= min_coverage, total_alt / cov, 0)
    
    allele_freq[i, ] <- af
  }
  
  return(allele_freq)
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

call_somatic_variants <- function(mgatk_data,
                                  seurat_object,
                                  baseline_cells = NULL,
                                  baseline_column = NULL,
                                  baseline_values = NULL,
                                  test_cells = NULL,
                                  test_column = NULL, 
                                  test_values = NULL,
                                  min_coverage = 10,
                                  min_alt_reads = 3,
                                  min_strand_reads = 2,
                                  min_cells_with_mutation = 3,
                                  max_baseline_af = 0.01,
                                  min_somatic_af = 0.05,
                                  strand_bias_max = 0.9,
                                  fisher_pval_cutoff = 0.01,
                                  p_adjust_method = "fdr") {
  
  if (is.null(baseline_cells) && !is.null(baseline_column)) {
    baseline_cells <- seurat_object@meta.data %>%
      tibble::rownames_to_column("barcode") %>%
      filter(get(baseline_column) %in% baseline_values) %>%
      pull(barcode)
  }
  
  if (is.null(test_cells) && !is.null(test_column)) {
    test_cells <- seurat_object@meta.data %>%
      tibble::rownames_to_column("barcode") %>%
      filter(get(test_column) %in% test_values) %>%
      pull(barcode)
  }
  
  if (is.null(baseline_cells) || is.null(test_cells)) {
    stop("Must specify baseline_cells and test_cells or provide column filters")
  }
  
  # Match barcodes with mgatk data
  baseline_barcodes <- intersect(baseline_cells, mgatk_data$barcodes)
  test_barcodes <- intersect(test_cells, mgatk_data$barcodes)
  
  
  if (length(baseline_barcodes) == 0 || length(test_barcodes) == 0) {
    stop("No matching cells found in mgatk data")
  }
  
  # Get cell indices
  baseline_idx <- which(mgatk_data$barcodes %in% baseline_barcodes)
  test_idx <- which(mgatk_data$barcodes %in% test_barcodes)
  
  positions <- mgatk_data$positions
  refallele <- mgatk_data$refallele
  bases <- c("A", "C", "G", "T")
  
  somatic_variants <- list()
  
  pb <- txtProgressBar(min = 0, max = length(bases), style = 3)
  
  for (base_idx in seq_along(bases)) {
    base <- bases[base_idx]
    
    fwd <- mgatk_data$counts[[paste0(base, "_fwd")]]
    rev <- mgatk_data$counts[[paste0(base, "_rev")]]
    cov <- mgatk_data$counts$coverage
    
    # Only test non-reference positions
    is_ref <- refallele == base
    alt_positions <- which(!is_ref)
    
    for (pos_idx in alt_positions) {
      pos <- positions[pos_idx]
      ref_base <- refallele[pos_idx]
      
      # Extract counts for this position
      fwd_counts <- fwd[, pos_idx]
      rev_counts <- rev[, pos_idx]
      cov_counts <- cov[, pos_idx]
      total_alt <- fwd_counts + rev_counts
      
      # Calculate allele frequencies
      af <- ifelse(cov_counts >= min_coverage, total_alt / cov_counts, NA)
      
      # Strand bias calculation
      strand_bias <- ifelse(total_alt > 0, 
                            pmax(fwd_counts, rev_counts) / total_alt, 0)
      
      # Quality filters per cell
      high_quality <- (
        cov_counts >= min_coverage &
          total_alt >= min_alt_reads &
          fwd_counts >= min_strand_reads &
          rev_counts >= min_strand_reads &
          strand_bias <= strand_bias_max &
          !is.na(af)
      )
      
      # Baseline statistics
      baseline_af <- af[baseline_idx]
      baseline_hq <- high_quality[baseline_idx]
      baseline_coverage <- cov_counts[baseline_idx]
      baseline_alt <- total_alt[baseline_idx]
      
      # Calculate baseline allele frequency (key for somatic calling)
      baseline_cells_covered <- sum(baseline_coverage >= min_coverage)
      if (baseline_cells_covered < 3) next
      
      baseline_mean_af <- mean(baseline_af[!is.na(baseline_af)])
      baseline_bulk_af <- sum(baseline_alt) / sum(baseline_coverage)
      baseline_cells_with_variant <- sum(baseline_hq & baseline_af > 0)
      
      # SOMATIC FILTER 1: Must be rare/absent in baseline
      if (baseline_bulk_af > max_baseline_af) next
      
      # Test population statistics  
      test_af <- af[test_idx]
      test_hq <- high_quality[test_idx]
      test_coverage <- cov_counts[test_idx]
      test_alt <- total_alt[test_idx]
      
      test_cells_covered <- sum(test_coverage >= min_coverage)
      if (test_cells_covered < min_cells_with_mutation) next
      
      test_mean_af <- mean(test_af[!is.na(test_af)])
      test_bulk_af <- sum(test_alt) / sum(test_coverage)
      test_cells_with_variant <- sum(test_hq & test_af > 0)
      
      # SOMATIC FILTER 2: Must be present in test population
      if (test_bulk_af < min_somatic_af) next
      if (test_cells_with_variant < min_cells_with_mutation) next
      
      # Statistical test: baseline vs test
      baseline_valid <- baseline_af[!is.na(baseline_af)]
      test_valid <- test_af[!is.na(test_af)]
      
      # Fisher's exact test on count data (more appropriate)
      fisher_pval <- NA
      if (length(baseline_valid) >= 3 && length(test_valid) >= 3) {
        cont_table <- matrix(c(
          sum(baseline_alt), sum(baseline_coverage) - sum(baseline_alt),
          sum(test_alt), sum(test_coverage) - sum(test_alt)
        ), nrow = 2, byrow = TRUE)
        
        tryCatch({
          fisher_result <- fisher.test(cont_table)
          fisher_pval <- fisher_result$p.value
        }, error = function(e) fisher_pval <<- NA)
      }
      
      if (is.na(fisher_pval) || fisher_pval > fisher_pval_cutoff) next
      
      # Calculate confidence intervals for effect size
      af_difference <- test_bulk_af - baseline_bulk_af
      odds_ratio <- (test_bulk_af / (1 - test_bulk_af)) / 
        (baseline_bulk_af / (1 - baseline_bulk_af) + 1e-10)
      
      # Allele frequency in cells with high-quality calls
      hq_baseline_af <- mean(baseline_af[baseline_hq], na.rm = TRUE)
      hq_test_af <- mean(test_af[test_hq], na.rm = TRUE)
      
      # Overall strand correlation
      all_fwd <- c(fwd_counts[baseline_idx], fwd_counts[test_idx])
      all_rev <- c(rev_counts[baseline_idx], rev_counts[test_idx])
      has_variant_idx <- (all_fwd + all_rev) > 0
      
      strand_correlation <- if (sum(has_variant_idx) > 1) {
        cor(all_fwd[has_variant_idx], all_rev[has_variant_idx], 
            use = "complete.obs")
      } else 0
      if (is.na(strand_correlation)) strand_correlation <- 0
      
      somatic_variants[[length(somatic_variants) + 1]] <- tibble(
        chromosome = "chrM",  # Mitochondrial
        position = pos,
        ref_allele = ref_base,
        alt_allele = base,
        variant_id = paste0("chrM:", pos, ":", ref_base, ">", base),
        
        # Baseline (normal) population
        baseline_cells = length(baseline_barcodes),
        baseline_coverage = sum(baseline_coverage),
        baseline_alt_reads = sum(baseline_alt),
        baseline_af = baseline_bulk_af,
        baseline_mean_af = baseline_mean_af,
        baseline_cells_with_variant = baseline_cells_with_variant,
        
        # Test (tumor/case) population  
        test_cells = length(test_barcodes),
        test_coverage = sum(test_coverage),
        test_alt_reads = sum(test_alt),
        test_af = test_bulk_af,
        test_mean_af = test_mean_af,
        test_cells_with_variant = test_cells_with_variant,
        
        # Somatic calling metrics
        af_difference = af_difference,
        fold_change = test_bulk_af / (baseline_bulk_af + 1e-10),
        odds_ratio = odds_ratio,
        fisher_pval = fisher_pval,
        
        # Quality metrics
        strand_correlation = strand_correlation,
        hq_baseline_af = hq_baseline_af,
        hq_test_af = hq_test_af
      )
    }
    
    setTxtProgressBar(pb, base_idx)
  }
  close(pb)
  
  if (length(somatic_variants) == 0) {
    message("No somatic variants detected")
    return(tibble())
  }
  
  # Combine results and apply multiple testing correction
  results <- bind_rows(somatic_variants) %>%
    mutate(
      fisher_pval_adj = p.adjust(fisher_pval, method = p_adjust_method),
      is_somatic = fisher_pval_adj < 0.05,
      somatic_score = -log10(fisher_pval_adj) * af_difference
    ) %>%
    arrange(fisher_pval)
  
  return(results)
}
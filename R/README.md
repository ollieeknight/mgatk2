# R helpers for mgatk2 HDF5 output

Requires hdf5r, Matrix, dplyr, and tibble; the example script also uses ggplot2.

```r
source("R/mgatk2_functions.R")

data <- read_mgatk_hdf5("path/to/output_directory")
data <- filter_cells_by_coverage(data, min_mean_depth = 10, min_coverage_breadth = 0.9)
data <- recompute_reference_alleles(data)

variants <- identify_variants(data, min_cells = 5, min_strand_cor = 0.65, min_vmr = 0.01)
allele_freq <- calculate_allele_freq(data, variants)
```

`R/mgatk2_qc_plots.R` runs the same steps with QC plots.

`read_mgatk_hdf5` returns a list of:

- `mats`: sparse cells x positions matrices `A_fwd` ... `T_rev`, `tn5_cuts_fwd`, `tn5_cuts_rev`, and `coverage`
- `cells`: one row per cell with `barcode`, `mean_depth`, `median_depth`, `max_depth`, `coverage_breadth`, `total_bases`, and any `barcode_metadata` columns
- `refallele`: one base per position, `N` where nothing was observed

## Functions

| Function | Notes |
| --- | --- |
| `read_mgatk_hdf5(dir)` | Loads one run. |
| `rbind_mgatk_data(...)` | Combines runs on the same reference. |
| `subset_cells(data, keep)` | Keeps cells by logical or index vector. |
| `subset_mgatk_barcodes(data, barcodes)` | Keeps the named barcodes. |
| `filter_cells_by_coverage(data, min_mean_depth, min_coverage_breadth)` | Thresholds differ by assay, so both are required. |
| `recompute_reference_alleles(data)` | The reference is inferred from aggregate counts; refresh it after subsetting. |
| `position_stats(data)` | Per-position depth, dropout, strand bias, and Tn5 cuts. |
| `identify_variants(data, ...)` | mgatk/Signac statistics: mean, VMR, strand correlation, cells detected. |
| `calculate_allele_freq(data, variants)` | Sparse variants x cells allele frequencies. |

`identify_variants` follows Signac: VMR is the variance of per-cell allele
frequency over the bulk frequency plus `1e-11`; `stabilise_variance = TRUE`
(the default) holds cells below `low_coverage_threshold` at the bulk frequency;
strand correlation is taken over cells carrying the allele on either strand.
For RNA, where strand carries no signal, pass `min_strand_cor = -1`.

## HDF5 layout

Both files live in `output/` and store matrices as positions x cells, chunked
128 cells wide; hdf5r reads them transposed, as cells x positions.

- `counts.h5`: `A_fwd`, `A_rev`, `C_fwd`, `C_rev`, `G_fwd`, `G_rev`, `T_fwd`,
  `T_rev`, `tn5_cuts_fwd`, `tn5_cuts_rev` (`uint16`), and `barcode`.
  Attributes `n_cells`, `n_positions`, `mito_chr`.
- `metadata.h5`: `coverage` (`uint16`), per-cell `mean_depth`, `median_depth`,
  `genome_coverage`, `total_bases` (`float32`) and `max_depth` (`uint16`),
  `reference` (one byte per position), and an optional `barcode_metadata/`
  group. Attributes `mito_chr`, `mito_length`.

`genome_coverage` is a fraction, matching `coverage_breadth` in
`qc/cell_stats.csv`.

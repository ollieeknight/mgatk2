# R analysis scripts for mgatk2 HDF5 output

## Usage

```r
source("R/mgatk2_functions.R")
library(hdf5r)
library(dplyr)

mgatk_data <- read_mgatk_hdf5("path/to/output_directory")

variants <- identify_variants(mgatk_data, min_cells = 5, min_strand_cor = 0.65, min_vmr = 0.01)
allele_freq <- calculate_allele_freq(mgatk_data, variants)

close_mgatk_hdf5(mgatk_data)
```

`read_mgatk_hdf5` leaves the count matrices on disk and returns:

- `barcodes`, `positions` (`1:16569`), `refallele` (one base per position, `N` where the pipeline never saw coverage)
- `mean_depth`, `median_depth`, `max_depth`, `genome_coverage` (fraction of positions covered), `total_bases` — one value per cell, named by barcode
- `barcode_metadata` — tibble of per-cell metrics if the run had a barcode metadata file, otherwise `NULL`
- `counts`, `metadata` — the open `H5File` handles
- `cells` — indices into the files, so `subset_mgatk_barcodes` and `filter_cells_by_coverage` are free

Call `close_mgatk_hdf5` when done.

Every function that touches a matrix streams it in blocks of 512 cells, which is how the files are chunked. Peak memory is a few hundred MB regardless of run size. Set `mgatk_data$block_size` to change it.

## Functions

| Function | Notes |
| --- | --- |
| `read_mgatk_hdf5(dir)` | Opens the pair of HDF5 files. |
| `close_mgatk_hdf5(data)` | Closes them. |
| `subset_mgatk_barcodes(data, barcodes)` | Restricts to matching barcodes; no data is read. |
| `filter_cells_by_coverage(data, min_mean_coverage, min_coverage_breadth)` | Thresholds are required — the sensible values differ by assay. |
| `calculate_cell_coverage_stats(data, deep = FALSE)` | Uses the pipeline's own per-cell numbers; `deep = TRUE` adds the 10x/50x counts and CV, which need a full pass. |
| `calculate_position_coverage_stats(data)` | Per-position depth summary. |
| `calculate_strand_coverage_stats(data)` | Forward/reverse balance per position. |
| `calculate_transposition_stats(data)` | Tn5 cut sites per position. |
| `identify_variants(data, ...)` | One streaming pass; returns the usual mgatk statistics. |
| `calculate_allele_freq(data, variants)` | Dense variants x cells matrix; `as(x, "sparseMatrix")` if you need it sparse. |

`identify_variants` follows Signac's conventions: the strand correlation uses cells carrying the alt allele on either strand, and VMR is the variance of the per-cell allele frequency over the bulk frequency plus a `1e-11` pseudocount. `stabilise_variance = TRUE` holds cells below `low_coverage_threshold` at the bulk frequency.

## Tests

```bash
Rscript R/test_mgatk2_functions.R
```

Builds a small HDF5 pair with a planted variant and asserts it comes back with the right frequency, strand correlation and VMR. Run from the repository root.

## HDF5 file structure

Both files live in `output/` and store matrices as (16569 positions x cells); hdf5r reads them back transposed, as (cells x 16569 positions).

### counts.h5

`A_fwd`, `A_rev`, `C_fwd`, `C_rev`, `G_fwd`, `G_rev`, `T_fwd`, `T_rev`, `tn5_cuts_fwd`, `tn5_cuts_rev` (all `uint16`), plus `barcode` (string array, note the singular name). Attributes: `n_cells`, `n_positions`, `mito_chr`.

### metadata.h5

`coverage` (`uint16`), the per-cell vectors `mean_depth`, `median_depth`, `genome_coverage`, `total_bases` (`float32`) and `max_depth` (`uint16`), `reference` (`S1`, one byte per position), and the optional `barcode_metadata/` group with one dataset per metadata column. Attributes: `mito_chr`, `mito_length`.

Matrices are chunked at 128 cells x all positions, so reading blocks of cells is cheap and reading single positions is not.

`genome_coverage` is a fraction, matching the `coverage_breadth` column of `qc/cell_stats.csv`. Files written before this change (including `scrna_realdata_output/`) store it as a percentage; divide by 100 before passing thresholds to `filter_cells_by_coverage`.

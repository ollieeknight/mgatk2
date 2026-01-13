# R analysis scripts for mgatk2 HDF5 output

## Files

- **mgatk_functions.R** - Data loading and variant calling functions
- **mgatk_qc_plots.R** - Line-by-line QC plotting script

## Usage

```r
source("R/mgatk_functions.R")
library(hdf5r)
library(tidyverse)

# Load HDF5 data
mgatk_data <- read_mgatk_hdf5("path/to/output_directory")

# Complete data structure:
# mgatk_data$barcodes: cell barcodes (character vector)
# mgatk_data$positions: 1:16569 (mtDNA positions)
# mgatk_data$refallele: reference allele at each position (character vector)
# mgatk_data$mean_coverage: mean coverage per position per cell (named numeric vector)
# mgatk_data$total_bases: total aligned bases per cell (named numeric vector)
# mgatk_data$depth: DEPRECATED - use mean_coverage instead
# mgatk_data$barcode_metadata: tibble with singlecell.csv metrics (if available)
# mgatk_data$counts: list containing count matrices (cells × positions in R):
#   - A_fwd, A_rev: A nucleotide counts (forward/reverse strand)
#   - C_fwd, C_rev: C nucleotide counts (forward/reverse strand)
#   - G_fwd, G_rev: G nucleotide counts (forward/reverse strand)
#   - T_fwd, T_rev: T nucleotide counts (forward/reverse strand)
#   - tn5_cuts_fwd: forward strand Tn5 cut sites (transposition events)
#   - tn5_cuts_rev: reverse strand Tn5 cut sites (transposition events)
#   - coverage: per-position coverage (cells × positions in R)

# Identify variants (slow for large datasets)
variants <- identify_variants(mgatk_data, min_cells = 5, min_strand_cor = 0.65, min_vmr = 0.01)

# Calculate allele frequencies
allele_freq <- calculate_allele_freq(mgatk_data, variants$variant)
```

## HDF5 File Structure

The mgatk2 output consists of two HDF5 files in the `output/` directory:

### counts.h5
Contains nucleotide count matrices stored as (16569 positions × cells) in Python/HDF5 format.
When read in R with hdf5r, matrices are transposed to (cells × 16569 positions):
- `A_fwd`, `A_rev`: A nucleotide counts (uint16)
- `C_fwd`, `C_rev`: C nucleotide counts (uint16)  
- `G_fwd`, `G_rev`: G nucleotide counts (uint16)
- `T_fwd`, `T_rev`: T nucleotide counts (uint16)
- `tn5_cuts_fwd`, `tn5_cuts_rev`: Tn5 transposase cut sites (uint16)
- `barcodes`: Cell barcodes (string array)

### metadata.h5
Contains metadata and summary information:
- `coverage`: Per-position coverage matrix (uint16, stored as positions × cells, read as cells × positions in R)
- `mean_coverage`: Mean coverage per position per cell (float32, length = number of cells)
- `total_bases`: Total aligned bases per cell (float32, length = number of cells)
- `reference`: Reference sequence (S1 bytes, length 16569, one per position)
- `barcode_metadata/`: Group containing per-cell QC metrics from singlecell.csv (if available):
  - `barcode`: Cell barcodes
  - `total_fragments`: Total fragments per cell
  - `mitochondrial_fragments`: Mitochondrial fragments per cell
  - `passed_filters`: Fragments passing QC filters
  - `peak_region_fragments`: Fragments in peaks
  - `TSS_fragments`: Fragments at TSS regions
  - Additional 10x ATAC-seq QC metrics as defined in singlecell.csv
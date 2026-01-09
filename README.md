# mgatk2: Mitochondrial Genome Analysis Toolkit v2

**mgatk2** is a reimplementation of the [original mgatk](https://github.com/caleblareau/mgatk) toolkit by [Caleb Lareau](https://github.com/caleblareau), optimised for processing mitochondrial DNA from single-cell ATAC-seq and other single-cell sequencing data. Tested with datasets of 200M chrM reads across 10k cells, this pipeline works well with datasets of all sizes.

## Key Improvements

- **Direct FASTA file hardmasking**: Critical for accurately mapping mitochondrial reads prior to running cellranger-atac
- **Pure Python implementation**: Single-pass processing without Snakemake or Java dependencies
- **HDF5 output format**: Fast binary storage with compression for efficient data access and analysis in R/Python
- **Stringent deduplication**: Fragment-length aware deduplication (default)
- **Tn5 cut site tracking**: Tracking of true transposition events
- **Optimised performance**: Parallel processing with efficient memory management
- **HTML QC reports**: Automatic generation of quality control visualizations (when singlecell.csv is available)

## Installation

```bash
# Create conda environment with Python 3.12 (recommended)
conda create -y -n mgatk2 python=3.12
conda activate mgatk2
pip install git+https://github.com/ollieeknight/mgatk2.git
```

## Quick Start

### Running with 10x Genomics data

**For Signac/Seurat compatibility** (text output format):
```bash
cd path/to/cellranger/output/  # directory containing 'outs'
mgatk2 tenx
```
This uses the original mgatk defaults (no quality filtering, alignment-only deduplication) and produces text output compatible with Signac's `ReadMGATK()`, `IdentifyVariants()` and `AlleleFreq()` functions.

**For enhanced QC and fast HDF5 output**:
```bash
cd path/to/cellranger/output/  # directory containing 'outs'
mgatk2 run
```
This enables stringent quality filtering, fragment-length aware deduplication, and generates HDF5 output for fast analysis with the included R functions (`R/mgatk2_functions.R`).

### Testing your configuration

```bash
mgatk2 test --dry-run  # Preview configuration without running
mgatk2 test            # Run on subsampled data for quick validation
```

## Output Files

### HDF5 format (default for `mgatk2 run`)

```
output/
├── counts.h5          # Nucleotide counts, Tn5 cut sites (16569 positions × cells)
├── metadata.h5        # Coverage, depth, reference, barcode metadata
└── report.html        # QC plots (generated if singlecell.csv is present)
```

### Text format (for `mgatk2 tenx` or `--format txt`)

```
output/
├── A.txt.gz, C.txt.gz, G.txt.gz, T.txt.gz  # Nucleotide count matrices
├── coverage.txt.gz                          # Coverage matrix
├── depthTable.txt                           # Depth per cell
├── chrM_refAllele.txt                       # Reference alleles
└── ...
```

## Command Reference

### Common options (all commands)

```bash
-i, --input PATH           Input BAM file or 10x directory, autodetected if 10X scATAC-seq run
-o, --output PATH          Output directory (default mgatk2)
-bt, --barcode-tag TEXT    BAM tag for cell barcodes (default CB)
-b, --barcodes PATH        Barcode file, , autodetected from outs/singlecell.csv if 10X scATAC-seq run
-g, --genome TEXT          Mitochondrial chromosome name (default chrM)
-t, --threads INTEGER      Number of parallel threads (default all)
--dry-run                  Show configuration without running
```

### Quality control options (`mgatk2 run` only)

```bash
-q, --quality INTEGER      Minimum base quality (Phred, default 20)
--mapq INTEGER             Minimum mapping quality (default 30)
-c, --min-reads INTEGER    Minimum reads per cell (default 10)
-s, --max-strand-bias      Maximum strand bias 0-1 (default 1.0)
-e, --min-distance-from-end  Min distance from read ends (default 5)
-d, --deduplication        Strategy: alignment_and_fragment_length (default),
                          alignment_start, or none
-f, --format              Output format: hdf5 (default) or txt
```

## Deduplication Methods

**mgatk2** provides three deduplication strategies:

1. **`alignment_and_fragment_length`** (default for `mgatk2 run`): Marks reads as duplicates if they share the same alignment position, strand, AND fragment length. Most stringent method, best for paired-end data.

2. **`alignment_start`** (default for `mgatk2 tenx`): Marks reads as duplicates if they share the same alignment position and strand only. Similar to Picard MarkDuplicates, matches original mgatk behavior.

3. **`none`**: No deduplication.

## License

MIT License - see [LICENSE](LICENSE) for details.

mgatk2 is derived from the original [mgatk](https://github.com/caleblareau/mgatk) by Caleb Lareau (2020).

## Citation

If you use mgatk2 in your research, please cite the original mgatk paper that established the methodology:

> Lareau CA, Ludwig LS, Muus C, et al. Massively parallel single-cell mitochondrial DNA genotyping and chromatin profiling. *Nature Biotechnology* 39, 451–461 (2021). https://doi.org/10.1038/s41587-020-0645-6
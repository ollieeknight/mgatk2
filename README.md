# mgatk2: Mitochondrial Genome Analysis Toolkit v2

**mgatk2** is a reimplementation of the [original mgatk](https://github.com/caleblareau/mgatk) toolkit by [Caleb Lareau](https://github.com/caleblareau), optimised for processing mitochondrial DNA from single-cell ATAC-seq and other single-cell sequencing data. Tested with datasets of 200M chrM reads across 10k cells, this pipeline works well with datasets of all sizes.

## Key Improvements

- **FASTA file hardmasking**: Mask regions of the genome highly similar to mitochondrial DNA, improving mapping quality
- **Python implementation**: Removed Snakemake and Java dependencies
- **HDF5 output format**: Fast binary storage with compression for efficient data access and analysis in R/Python
- **Stringent deduplication**: Fragment-length aware deduplication alongside position and strand as in original
- **Optimised performance**: Parallel processing with efficient memory management
- **HTML QC reports**: Automatic generation of quality control visualisations

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
cd path/to/cellranger-atac/output/  # directory containing 'outs'
mgatk2 tenx
```
This uses the original mgatk defaults (no quality filtering, alignment-only deduplication) and produces text output compatible with Signac's `ReadMGATK()`, `IdentifyVariants()` and `AlleleFreq()` functions.

**For enhanced QC and fast HDF5 output**:
```bash
cd path/to/cellranger-atac/output/  # directory containing 'outs'
mgatk2 run
```
This enables stringent quality filtering, fragment-length aware deduplication, and generates HDF5 output for fast analysis with the included R functions (`R/mgatk2_functions.R`).

## Output Files

### HDF5 format (default for `mgatk2 run`)

```
tree tests/run_output/
├── mgatk2_report.html     # Quality control report
├── output
│   ├── counts.h5          # Nucleotide counts, Tn5 cut sites (16569 positions × cells)
│   └── metadata.h5        # Coverage, depth, reference, barcode metadata
├── output.log             # Complete log of run command
└── qc
    ├── cell_stats.csv     # Cell by cell statistics in a .csv format
    └── summary.txt        # Quick output summary
```

### Text format (for `mgatk2 tenx` or `--format txt`)

```
tree tests/tenx_output
├── output
│   ├── A.txt.gz              # Nucleotide count matrix
│   ├── C.txt.gz              # Nucleotide count matrix
│   ├── chrM_refAllele.txt    # Reference alleles
│   ├── coverage.txt.gz       # Coverage matrix
│   ├── depthTable.tx         # Depth per cell
│   ├── G.txt.gz              # Nucleotide count matrix
│   └── T.txt.gz              # Nucleotide count matrix
├── output.log                # Complete log of run command
└── qc
    ├── cell_stats.csv        # Cell by cell statistics in a .csv format
    └── summary.txt           # Quick output summary

```

## Command Reference

```bash
mgatk2 run --help
Usage: mgatk2 run [OPTIONS]

  Run mgatk2 with optimised defaults

Options:
  -i, --input PATH                      Input BAM file or 10x outs/ directory
                                        [default: current directory, auto-detects possorted_bam.bam]
  -g, --genome TEXT                     Mitochondrial chromosome name (e.g chrM, MT, or M)  
                                        [default: chrM]
  -b, --barcodes PATH                   Barcode file (singlecell.csv, barcodes.tsv/csv, or auto-detect from BAM)
  -bt, --barcode-tag TEXT               BAM tag for cell barcode  
                                        [default: CB]
  --min-barcode-reads INTEGER           Minimum reads per barcode when auto-detecting from BAM  
                                        [default: 10]
  -o, --output PATH                     Output directory for analysis results  
                                        [default: mgatk2/]
  -t, --threads INTEGER                 Number of threads for parallel processing 
                                        [default: 16]
  -v, --verbose                         Enable verbose logging 
                                        [default: on]
  --batch-size INTEGER                  Number of cells to process per batch
                                        [default: 250]
  -m, --memory FLOAT                    Maximum memory usage in GB 
                                        [default: 128]
  -q, --quality INTEGER                 Minimum base quality (Phred score) 
                                        [default: 20]
  --mapq INTEGER                        Minimum alignment/mapping quality  
                                        [default: 30]
  -c, --min-reads INTEGER               Minimum deduplicated reads per cell to include in analysis  
                                        [default: 1]
  -s, --max-strand-bias FLOAT           Maximum strand bias (0-1)  
                                        [default: 1.0]
  -e, --min-distance-from-end INTEGER   Minimum distance from read ends (bp) 
                                        [default: 5]
  -d, --deduplication                   Deduplication strategy 
                                        [alignment_and_fragment_length|alignment_start|none]
                                        [default: alignment_and_fragment_length]
  -f, --format [txt|hdf5]               Output format: txt (text files) or hdf5 (fast binary)
                                        [default: hdf5]
  --dry-run                             Show configuration and exit without processing (no output files created)
  --help                                Show this message and exit.
```

### `mgatk2 call` - Bulk Analysis

```bash
mgatk2 call --help
Usage: mgatk2 call [OPTIONS]

  Bulk analysis of BAM files (Smart-seq style, one BAM per cell)

Options:
  -i, --input PATH                      Input BAM file for bulk analysis (one BAM per cell)  [required]
  -g, --genome TEXT                     Mitochondrial chromosome name (e.g chrM, MT, or M)  
                                        [default: chrM]
  -o, --output PATH                     Output directory for analysis results  
                                        [default: mgatk2/]
  -t, --threads INTEGER                 Number of threads for parallel processing 
                                        [default: 16]
  -v, --verbose                         Enable verbose logging 
                                        [default: on]
  -m, --memory FLOAT                    Maximum memory usage in GB 
                                        [default: 128]
  -q, --quality INTEGER                 Minimum base quality (Phred score) 
                                        [default: 20]
  --mapq INTEGER                        Minimum alignment/mapping quality  
                                        [default: 30]
  -s, --max-strand-bias FLOAT           Maximum strand bias (0-1)  
                                        [default: 1.0]
  -e, --min-distance-from-end INTEGER   Minimum distance from read ends (bp) 
                                        [default: 5]
  -d, --deduplication                   Deduplication strategy 
                                        [alignment_and_fragment_length|alignment_start|none]
                                        [default: alignment_and_fragment_length]
  -f, --format [txt|hdf5]               Output format: txt (text files) or hdf5 (fast binary)
                                        [default: hdf5]
  --dry-run                             Show configuration and exit without processing
  --help                                Show this message and exit.
```

### `mgatk2 hardmask-fasta` - FASTA Masking

```bash
mgatk2 hardmask-fasta --help
Usage: mgatk2 hardmask-fasta [OPTIONS]

  Hard-mask reference genome FASTA with blacklists

Options:
  -i, --input-fasta PATH                Input reference genome FASTA file  [required]
  -o, --output-fasta PATH               Output hard-masked FASTA file  [required]
  -g, --genome TEXT                     Genome build (hg38, hg19, GRCh38, GRCh37, mm10, mm9, GRCm38, GRCm37)  [required]
  --mt-chrom TEXT                       Mitochondrial chromosome name (auto-detected from genome if not provided)
  --mask-numts / --no-mask-numts        Mask NUMT regions in nuclear chromosomes (recommended)  
                                        [default: mask-numts]
  -v, --verbose                         Enable verbose logging  
                                        [default: off]
  --help                                Show this message and exit.
```

## Deduplication Methods

**mgatk2** provides three deduplication strategies:

1. **`alignment_and_fragment_length`** (default for `mgatk2 run`): Marks reads as duplicates if they share the same alignment position, strand, AND fragment length. Most stringent method and the best for paired-end data.

2. **`alignment_start`** (default for `mgatk2 tenx`): Marks reads as duplicates if they share the same alignment position and strand only. This mode is similar to `Picard MarkDuplicates` and matches original mgatk behavior.

3. **`none`**: No deduplication.

## License

MIT License - see [LICENSE](LICENSE) for details.

mgatk2 is derived from the original [mgatk](https://github.com/caleblareau/mgatk) by Caleb Lareau (2020).

## Citation

If you use mgatk2 in your research, please cite the original mgatk paper that established the methodology:

> Lareau CA, Ludwig LS, Muus C, et al. Massively parallel single-cell mitochondrial DNA genotyping and chromatin profiling. *Nature Biotechnology* 39, 451–461 (2021). https://doi.org/10.1038/s41587-020-0645-6
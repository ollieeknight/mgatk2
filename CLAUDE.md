# mgatk2 — CLAUDE.md

Mitochondrial genome analysis toolkit for single-cell data. Rewrite of Caleb Lareau's mgatk with a focus on performance and HDF5-native output.

## Project structure

```
src/
  cli/           # Click CLI — base.py, options.py, commands/{run,tenx,call,mask}.py
  core/          # config.py (PipelineConfig dataclasses), pipeline.py (MtDNAPipeline), exceptions.py
  processing/    # pileup.py (PileupGenerator), processors.py (CellProcessor), readers.py (BAMReader)
  file_io/       # writers.py (IncrementalHDF5Writer, IncrementalTextWriter), formats.py, main.py
  analysis/      # qc.py, report.py (HTML report generation)
  utils/         # genome_utils.py, masking.py, utils.py
  data/blacklists/  # hg38/hg19/mm10/mm9 NUMTs BED files
R/
  mgatk2_functions.R   # read_mgatk_hdf5, subset, filter, coverage stats, variant calling
  mgatk2_qc_plots.R    # example QC plot workflow
tests/
  run_output/          # HDF5 output from mgatk2 run
  run_hdf5_output/     # same, explicit HDF5 mode
  run_txt_output/      # txt format output
  tenx_output/         # 10x/Signac-compatible txt output
```

## CLI commands

| Command | Use case | Key defaults |
|---------|----------|-------------|
| `mgatk2 run` | scATAC/scRNA BAM with cell barcodes | `--format hdf5`, dedup `alignment_and_fragment_length` |
| `mgatk2 tenx` | 10x Genomics outs/ directory | `--format txt`, dedup `alignment_start` (Signac-compatible) |
| `mgatk2 call` | Bulk BAM, one file per cell | `--format hdf5` |
| `mgatk2 mask` | Generate position mask BED | |

Common options across all: `--input/-i`, `--genome/-g` (chrM name), `--output/-o`, `--threads/-t`, `--format/-f`, `--deduplication/-d`, `--quality/-q`, `--mapq`, `--min-reads/-c`, `--max-strand-bias/-s`, `--min-distance-from-end/-e`, `--dry-run`.

## mgatk2 wes — somatic mitochondrial variant calling from WES CRAMs

New command added in v1.1. Takes matched tumour + normal CRAM files (e.g. from CLONK-WES) and
calls somatic mitochondrial variants without GATK, using the same pysam pileup engine as the
single-cell commands.

### Usage

```bash
mgatk2 wes \
  --tumour K010_NKG2C.cram \
  --normal K010_CD14.cram \
  --reference GRCh38.fa \
  --output mgatk2_wes/K010_NKG2C/ \
  --sample-name K010_NKG2C \
  --genome chrM \
  --min-af 0.005 \
  --max-normal-af 0.01 \
  --min-tn-ratio 3.0 \
  --min-depth 10 \
  --min-normal-depth 5 \
  --mapq 20
```

### Key design decisions

- **No GATK** — direct pysam pileup; no shifted-reference trick needed (only required for GATK's HMM model, not pileup)
- **CRAM support** — `pysam.AlignmentFile(path, "rc", reference_filename=ref)` for CRAM decoding; BAM also accepted
- **NuMT defence** — primary filter is `--mapq 20` (bwa-mem2 reduces MAPQ for reads ambiguously mapping to nuclear NuMTs vs chrM)
- **Bundled NuMT BEDs** (`src/data/blacklists/`) are **nuclear-side** positions used for reference masking, **not chrM-side** — they contain no chrM rows; blacklist default is `"none"`. Supply `--custom-blacklist` for chrM-side exclusions.
- **Reference alleles** inferred from majority base across both pileups (no separate FASTA required)
- **Deduplication disabled** — CLONK-WES CRAMs are already fgbio consensus-deduplicated
- **All candidate sites written** (not just PASS) — downstream R analysis can refilter

### Output

Single TSV: `<output>/<sample-name>.mito_somatic.tsv`

Columns: `CHROM POS REF ALT NORMAL_DEPTH TUMOUR_DEPTH NORMAL_REF_COUNT NORMAL_ALT_COUNT NORMAL_FWD NORMAL_REV TUMOUR_REF_COUNT TUMOUR_ALT_COUNT TUMOUR_FWD TUMOUR_REV NORMAL_VAF TUMOUR_VAF TN_RATIO STRAND_BIAS FILTER`

`FILTER` = `PASS` or pipe-delimited reasons:
`LOW_TUMOUR_DEPTH | LOW_NORMAL_DEPTH | LOW_ALT_READS | LOW_VAF | GERMLINE | LOW_TN_RATIO | STRAND_BIAS | BLACKLIST`

Log: `<output>/wes.log`

### New source files

| File | Role |
|------|------|
| `src/core/config.py:WesConfig` | Dataclass with all WES thresholds; `to_pipeline_config()` for PileupGenerator |
| `src/data/blacklists/__init__.py:load_blacklist_positions()` | Loads bundled or custom chrM-side BED → `set[int]` of 1-based positions |
| `src/processing/wes_pileup.py` | `CRAMReader`, `infer_reference_alleles`, `call_somatic_variants`, `write_variants_tsv`, `run_wes_pipeline` |
| `src/cli/options.py:wes_options()` | Click option decorator for `wes` command |
| `src/cli/commands/wes.py` | Click `wes` command |

### Filter defaults

| Filter | Default | Flag |
|--------|---------|------|
| Min tumour depth | 10 | `--min-depth` |
| Min normal depth | 5 | `--min-normal-depth` |
| Min tumour AF | 0.005 | `--min-af` |
| Max normal AF | 0.01 | `--max-normal-af` |
| Min T/N ratio | 3.0 | `--min-tn-ratio` |
| Min alt reads | 3 | `--min-alt-reads` |
| Max strand bias | 0.9 | `--max-strand-bias` |
| Min MAPQ | 20 | `--mapq` |

---

## Output formats

### HDF5 (default for `run` and `call`)

Two files in `<output>/output/`:

**`counts.h5`** — shape `(n_positions, n_barcodes)` in Python/C order; hdf5r reads as `(n_barcodes, n_positions)` in R (transposed automatically):
- `barcode` — array of barcode strings
- `{A,C,G,T}_fwd`, `{A,C,G,T}_rev` — base counts per strand, uint16
- `tn5_cuts_fwd`, `tn5_cuts_rev` — Tn5 insertion sites, uint16

**`metadata.h5`** — same shape convention:
- `coverage` — total read depth per position per cell, uint16 `(n_positions, n_barcodes)`
- `mean_depth`, `median_depth`, `max_depth`, `genome_coverage`, `total_bases` — per-barcode scalars
- `reference` — reference allele per position, `S1` byte array, length `n_positions`
- `barcode_metadata` (optional group) — from `singlecell.csv` (scATAC)

QC in `<output>/qc/`: `cell_stats.csv`, `summary.txt`.

### TXT (default for `tenx`, Signac-compatible)

Files in `<output>/output/`:
- `output.{A,C,G,T}.txt.gz` — `pos,barcode,fwd,rev` (no header, 1-based position)
- `output.coverage.txt.gz` — `pos,barcode,depth`
- `output.depthTable.txt` — `barcode\tdepth`
- `{chrM}_refAllele.txt` — `pos\tref` (has header)

## Key design decisions

- **Single-pass BAM processing**: reads collected by barcode then processed in parallel batches
- **HDF5 staged writes**: data written to local TMPDIR staging then moved to final location (network FS resilience)
- **Batch size**: `io_batch_size` defaults to `max(50, min(10% of barcodes, 1000))`; `worker_batch_size` defaults to `ncores`
- **Dedup modes**: `alignment_and_fragment_length` (start pos + fragment length), `alignment_start` (start only), `none`
- **Strand bias filter**: applied per-base per-position; positions exceeding `max_strand_bias` are zeroed (not dropped)
- **Tn5 cut sites**: recorded at read 5' end (fwd: `reference_start`, rev: `reference_start + read_length - 1`)
- **Reference allele**: computed post-hoc as most common base across all cells at each position
- **mito_length**: hardcoded 16569 (hg38 chrM); auto-detected from BAM header is not implemented — must match input

## R functions (mgatk2_functions.R)

All assume HDF5 format. Main entry: `read_mgatk_hdf5(dir)` returns a list with:
- `counts` — list of matrices (cells × positions): `A_fwd`, `A_rev`, ..., `tn5_cuts_fwd`, `tn5_cuts_rev`, `coverage`
- `mean_depth`, `median_depth`, `max_depth`, `genome_coverage`, `total_bases` — named vectors (by barcode)
- `refallele`, `positions`, `barcodes`, `barcode_metadata`

Downstream: `subset_mgatk_barcodes`, `filter_cells_by_coverage`, `calculate_cell_coverage_stats`, `calculate_position_coverage_stats`, `calculate_strand_coverage_stats`, `calculate_transposition_stats`, `recompute_reference_alleles`, `identify_variants`, `calculate_allele_freq`.

`identify_variants` returns tibble with: `position`, `nucleotide`, `variant` (e.g. `"302A>C"`), `mean` (bulk AF), `vmr`, `n_cells_detected`, `n_cells_conf_detected`, `strand_correlation`.

`calculate_allele_freq` returns a sparse `Matrix` (variants × cells).

## Dependencies

Python: `click>=8.0`, `pysam>=0.19`, `numpy>=1.20`, `tqdm>=4.50`, `h5py>=3.0`, `matplotlib>=3.3`

R: `hdf5r`, `tidyverse`, `Matrix` (core); `Signac`/`Seurat` for 10x integration; `pheatmap` for variant heatmaps

Container: `Dockerfile` at repo root; conda recipe at `conda-recipe/meta.yaml`

## Testing

Test data in `tests/outs/`: `possorted_bam.bam` + barcode files from a small 10x ATAC experiment. Expected outputs committed under `tests/run_output/`, `tests/run_hdf5_output/`, `tests/run_txt_output/`, `tests/tenx_output/`. CI via `.github/workflows/ci.yml`.

## Common gotchas

- `--genome` must match the chromosome name in the BAM header exactly (tries `chrM`, `MT`, `M`, `chrMT` in order)
- `mito_length` is hardcoded to 16569; non-human genomes need this adjusted in `PipelineConfig`
- `barcode_file` can be a `.csv` (singlecell.csv from 10x, triggers ATAC report), `.tsv`/`.txt` (plain list), or omitted (auto-extracted from BAM using `--barcode-tag`)
- HDF5 datasets are stored `(positions, cells)` in Python; hdf5r automatically transposes to `(cells, positions)` in R
- `tn5_cuts` are only meaningful for ATAC data; RNA/bulk will have zeros

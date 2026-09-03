# mgatk2: mitochondrial genome analysis toolkit

**mgatk2** is a reimplementation of the [original mgatk](https://github.com/caleblareau/mgatk)
toolkit by [Caleb Lareau](https://github.com/caleblareau), for mitochondrial DNA
from single-cell ATAC-seq and other single-cell sequencing data. It has been
tested on datasets of 200M chrM reads across 10k cells, and works comfortably at
both ends of that range.

Alongside the single-cell workflow it also hard-masks nuclear NUMT regions in a
reference FASTA, and compares a sample against an autologous normal to produce a
somatic-style mitochondrial SNV VCF.

## Key improvements

- **Pure Python**: single-pass processing, with no Snakemake and no Java
- **Stringent deduplication**: fragment-length aware by default
- **HDF5 output**: compressed binary matrices that load quickly in R or Python
- **Direct FASTA hard-masking**: essential for mapping mitochondrial reads
  honestly, before cellranger-atac ever runs
- **Tn5 cut site tracking**: one insertion recorded per retained read, so cut
  totals and read counts always agree
- **Bounded memory**: peak memory tracks the cell count and your budget, never
  the number of reads in the BAM
- **HTML QC reports**: written automatically for single-cell HDF5 runs
- **Paired tumour/normal calling**: a two-sample VCF that drops straight into
  `bcftools`, VEP, and the rest of the somatic tool chain

## Installation

mgatk2 needs Python 3.10 or newer.

```bash
# Conda environment (recommended)
conda create -y -n mgatk2 python=3.12
conda activate mgatk2
pip install git+https://github.com/ollieeknight/mgatk2.git
mgatk2 --help
```

Or with a plain virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
pip install git+https://github.com/ollieeknight/mgatk2.git
```

For development:

```bash
git clone https://github.com/ollieeknight/mgatk2.git
cd mgatk2
make setup
make
```

`make` formats, lints, runs the unit tests, and executes every single-cell
pipeline end to end against the committed fixtures.

## Quick start

### 10x Genomics data

**For Signac/Seurat compatibility** (text output):

```bash
cd path/to/cellranger/output/   # the directory containing 'outs'
mgatk2 tenx
```

This uses the original mgatk defaults, meaning no quality filtering and
alignment-only deduplication, and writes text output for Signac's `ReadMGATK()`,
`IdentifyVariants()`, and `AlleleFreq()`.

**For stringent QC and fast HDF5 output**:

```bash
cd path/to/cellranger/output/
mgatk2 run
```

This turns on quality filtering and fragment-length aware deduplication, writes
HDF5, and generates an HTML QC report. The R helpers in
[R/mgatk2_functions.R](R/mgatk2_functions.R) read that output; see
[R/README.md](R/README.md).

`--input` also accepts a BAM directly, or a 10x Multi / CITE-seq `outs/`
directory (`per_sample_outs/<sample>/count/sample_alignments.bam`). The matching
barcode file is found automatically in both cases. Where a Multi `outs/` holds
more than one demultiplexed sample, point `--input` at that sample's
`sample_alignments.bam`.

### Checking a configuration before committing to it

```bash
mgatk2 run --dry-run
```

Every command takes `--dry-run`. It prints the configuration, opens the input
alignments, confirms the mitochondrial contig and an index are both present, and
then exits without writing anything at all.

## Commands

| Command | Purpose |
|---|---|
| `mgatk2 run` | Single-cell BAM with barcode tags; filtered HDF5 defaults |
| `mgatk2 tenx` | 10x `outs/` input; text defaults matching original mgatk |
| `mgatk2 call` | A directory holding one bulk BAM per sample |
| `mgatk2 paired` | Tumour/normal mitochondrial SNV evidence from BAM or CRAM |
| `mgatk2 hardmask-fasta` | Hard-mask a reference FASTA with the bundled NUMT BEDs |

Every command accepts `-h` as well as `--help`, and
`mgatk2 COMMAND --help` is always the authority on the current defaults.

## Output files

### HDF5 (default for `run` and `call`)

```
output/
├── output/
│   ├── counts.h5        # A/C/G/T per strand, Tn5 cuts, barcodes (positions × cells)
│   └── metadata.h5      # coverage, per-cell depth summaries, inferred reference
├── qc/
│   ├── cell_stats.csv   # per-cell reads, fragments, depth, coverage breadth
│   ├── run_config.json  # the same record as JSON, for downstream scripts
│   └── summary.txt      # run summary and the parameters it ran with
├── mgatk2_report.html   # QC plots
└── output.log
```

Matrices are stored positions × cells; `hdf5r` reads them transposed, as cells ×
positions. Bulk runs write no HTML report, since it is a per-cell document and a
bulk sample is a single pseudo-cell.

### Text (default for `tenx`, or `--format txt`)

```
output/
├── output/
│   ├── output.A.txt.gz, output.C.txt.gz, output.G.txt.gz, output.T.txt.gz
│   ├── output.coverage.txt.gz
│   ├── output.depthTable.txt
│   └── chrM_refAllele.txt
├── qc/
└── output.log
```

The reference allele table is inferred from aggregate counts across all cells.

## Command reference

### Shared single-cell options (`run`, `tenx`, `call`)

```
-i, --input PATH             BAM file or 10x outs/ directory (a directory of BAMs for `call`)
-o, --output PATH            Output directory [default: mgatk2]
-g, --genome TEXT            Mitochondrial chromosome name [default: chrM]
-b, --barcodes PATH          Barcode file; auto-detected from a 10x outs/ directory
-bt, --barcode-tag TEXT      BAM tag holding the cell barcode [default: CB]
    --min-barcode-reads INT  Minimum reads per barcode when detecting them from the BAM [default: 10]
-t, --threads INT            Concurrent shard workers, or concurrent BAMs for `call` [default: auto]
-m, --memory FLOAT           Memory budget in GB shared by the workers [default: 128]
-q, --quality INT            Minimum base quality
    --mapq INT               Minimum mapping quality
-c, --min-reads INT          Minimum deduplicated reads per cell (floored at 1)
-s, --max-strand-bias FLOAT  Maximum per-base strand imbalance, 1.0 disables [default: 1.0]
-e, --min-distance-from-end INT   Ignore bases this close to either read end
    --nh-max INT             Maximum NH tag; 0 disables, 1 matches mgatk [default: 0]
    --nm-max INT             Maximum NM/nM tag; 0 disables, 4 matches mgatk [default: 0]
-d, --deduplication CHOICE   alignment_and_fragment_length | alignment_start | none
-f, --format CHOICE          hdf5 | txt
    --no-tn5                 Skip Tn5 cut site tracking, for RNA-seq or WGS
-v, --verbose                Verbose logging
    --dry-run                Validate and exit without processing
```

`call` has no barcode options: it treats every `*.bam` in the input directory as
one bulk sample, running one process per BAM so each keeps its own log file.
Its `--memory` is accepted for compatibility but can never bind, because a bulk
sample is a single cell and so never shards.

The defaults that differ between the presets:

| Option | `run` | `tenx` | `call` |
|---|---|---|---|
| `--quality` | 20 | 0 | 20 |
| `--mapq` | 30 | 0 | 30 |
| `--min-reads` | 1 | 0 | n/a |
| `--min-distance-from-end` | 5 | 0 | 5 |
| `--deduplication` | fragment length | alignment start | fragment length |
| `--format` | hdf5 | txt | hdf5 |

### `mgatk2 paired`

```
    --tumor FILE                 Tumour (query) BAM or CRAM [required]
    --normal FILE                Autologous normal (comparator) BAM or CRAM [required]
    --reference FILE             Indexed FASTA defining REF; also decodes CRAM [required]
-o, --output PATH                Output directory [required]
    --sample-name TEXT           Filename prefix for the VCF, index, and BED [required]
-g, --genome TEXT                Mitochondrial chromosome name [default: chrM]
-q, --quality INT                Minimum base quality [default: 20]
    --mapq INT                   Minimum mapping quality [default: 20]
-e, --min-distance-from-end INT  Discard observations this close to a read end [default: 5]
-s, --max-strand-bias FLOAT      Tumour alternate strand imbalance above this flags STRAND_BIAS [default: 0.9]
-d, --deduplication CHOICE       Coordinate fallback for unmarked input [default: fragment length]
    --min-tumor-depth INT        [default: 10]
    --min-normal-depth INT       [default: 5]
    --min-alt-observations INT   [default: 3]
    --min-tumor-af FLOAT         [default: 0.005]
    --max-normal-af FLOAT        [default: 0.01]
    --custom-blacklist FILE      BED of chrM positions to flag BLACKLIST
    --autosomal-median-depth FLOAT   Enables the POSSIBLE_NUMT filter
    --circular-edge-bases INT    Linear-reference edge width to flag [default: 500]
    --input-is-consensus         Declare UMI-consensus input (needs --deduplication none)
```

### `mgatk2 hardmask-fasta`

```
-i, --input-fasta FILE   Reference FASTA, plain or .gz [required]
-o, --output-fasta FILE  Output FASTA; gzipped when the name ends in .gz [required]
-g, --genome TEXT        hg38, hg19, GRCh38, GRCh37, mm10, mm9, GRCm38, GRCm37 [required]
```

```bash
mgatk2 hardmask-fasta -i GRCh38.fa -o GRCh38.numt_masked.fa.gz -g hg38
```

This masks nuclear NUMT regions only, so it takes no mitochondrial arguments.
Run it before building your alignment index: reads of nuclear mitochondrial
origin otherwise mismap onto chrM and quietly inflate heteroplasmy.

## Deduplication

Three strategies are available everywhere:

1. **`alignment_and_fragment_length`** (default for `run`, `call`, and `paired`):
   duplicates share alignment start, strand, *and* template length. The most
   stringent option, and the right one for paired-end data.
2. **`alignment_start`** (default for `tenx`): duplicates share alignment start
   and strand. Close to Picard MarkDuplicates, and matches original mgatk.
3. **`none`**: keep every otherwise eligible alignment. Use this for input that
   is already deduplicated or UMI-consensus collapsed.

Deduplication is applied per cell.

## How the single-cell scan works

The barcode list is split into contiguous shards. Each shard is handled by one
worker that opens the BAM itself, streams chrM once, and accumulates straight
into dense per-cell arrays. Reads are never held as Python objects, and a
finished shard travels back as NumPy blocks that land in HDF5 as one contiguous,
chunk-aligned column write.

`--threads` sets how many shard workers run at once, and `--memory` sets the
budget they share. Shard width follows from that budget, so a run with more
cells than the budget allows simply uses more, narrower shards rather than
growing its memory. Each shard decodes chrM independently, so decode work scales
with shard count while wall time does not: a deliberate trade for a single-pass
design with bounded memory and no inter-process read traffic.

Tn5 cut sites record one insertion per retained read, at the read's outermost
aligned reference base, so the cut total always equals the retained read count
in `qc/cell_stats.csv`. They are meaningful for ATAC data; use `--no-tn5` for
RNA-seq or WGS.

`--max-strand-bias` is the per-base imbalance `|forward - reverse| / total`, the
same metric `paired` uses. At the single-cell default of `1.0` the filter is off;
above the threshold, a base's counts are zeroed.

## Paired tumour/normal analysis

`paired` compares a sample against an autologous normal:

```bash
mgatk2 paired \
  --tumor tumour.bam \
  --normal normal.bam \
  --reference GRCh38.fa \
  --output results \
  --sample-name K010
```

The tumour/normal vocabulary is naming, not biology. Any two autologous samples
pair validly, and an NKG2C+ population against CD14+ from the same donor is the
same test. The output keeps somatic conventions so it slots into existing tools,
which is also why column names and VCF fields use American spelling (`tumor_af`,
`normal_dp`) while the prose stays British.

For UMI-consensus input, declare it and do not deduplicate a second time:

```bash
mgatk2 paired \
  --tumor tumour.consensus.cram \
  --normal normal.consensus.cram \
  --reference GRCh38.fa \
  --output results \
  --sample-name pair \
  --deduplication none \
  --input-is-consensus
```

The FASTA defines `REF`; observed major alleles never replace it. Overlapping
mates contribute at most one observation per position, with higher base quality
resolving a disagreement and equal-quality disagreements masked. Missing-quality,
duplicate, QC-failed, unmapped, secondary, and supplementary reads are excluded
and counted in QC.

### Outputs

| File | Contents |
|---|---|
| `{sample}.mt_variants.vcf.gz` | Two-sample SNV VCF, the complete variant record |
| `{sample}.mt_variants.vcf.gz.tbi` | Tabix index |
| `{sample}.mt_callable.bed.gz` | Positions reportable in both samples |

The callable BED earns its place because it is the one thing the VCF cannot
stand in for: a VCF holds variant sites, so callable territory would otherwise be
unrecoverable after the run. It is the denominator for any per-callable-base
mutation burden, and the same companion GATK CallableLoci and Strelka emit.

### The VCF is the record

Samples are ordered `NORMAL`, `TUMOR`.

- `FORMAT`: `GT:DP:AD:AF:SB:F1R2:F2R1:AFCI`, per sample.
- `INFO`: `DP`, `ERR`, `SEQP`, `EP`, `EQ`, `SP`, `OP`; the tumour-derived
  per-allele medians `MBQ`, `MMQ`, `MPOS` (`Number=R`, reference then alternate,
  as Mutect2 lays them out); and the alternate-versus-reference rank-sum
  z-scores and p-values `RSBQ`, `RSMQ`, `RSPOS`.
- `QUAL` is the phred-scaled, BH-adjusted enrichment q-value.

`AD` carries genuine reference and alternate counts rather than `DP - AC`, which
differ at a multiallelic position. Quality medians are tumour-derived only:
nothing filters on the normal-side values, and carrying them doubled every
record for no reader.

### QC and provenance live in the header

There is no sidecar JSON. The whole QC record — parameters, input paths and read
statistics, reference checksum, per-substitution error rates, depth quantiles,
candidate and callable-position counts, mgatk2 version, git commit, and the
command line — is one `##mgatk2_qc=` header line holding a single JSON object,
and it survives a `bcftools` round trip:

```bash
bcftools view -h sample.mt_variants.vcf.gz | sed -n 's/^##mgatk2_qc=//p' | jq
```

Read the variants with `VariantAnnotation::readVcf` in R, or flatten them:

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%DP\t%AD\t%AF]\n' \
  sample.mt_variants.vcf.gz
```

### How a candidate is judged

Two independent questions must both be answered before a variant passes:

1. **Is the tumour enriched over the normal?** A one-sided Fisher exact test,
   Benjamini-Hochberg adjusted across every observed allele (`EP`, `EQ`).
2. **Is the tumour alternate above sequencing noise?** A per-substitution error
   rate is learned from the normal, excluding any site whose normal allele
   fraction exceeds `--max-normal-af` because that is a plausible real allele,
   and the tumour alternate count is tested against it with a binomial (`ERR`,
   `SEQP`).

Both are needed. A Fisher test alone assumes the two samples share an error
rate, so a difference in depth can look significant on its own.

### Filters

| Filter | Fires when | Mutect2 analogue |
|---|---|---|
| `LOW_TUMOR_DEPTH` / `LOW_NORMAL_DEPTH` | Depth below `--min-tumor-depth` / `--min-normal-depth` | — |
| `LOW_ALT_OBSERVATIONS` | Tumour alternate count below `--min-alt-observations` | — |
| `LOW_TUMOR_AF` | Tumour AF below `--min-tumor-af` | `low_allele_frac` |
| `HIGH_NORMAL_AF` | Normal AF above `--max-normal-af` | `normal_artifact`, `germline` |
| `WEAK_EVIDENCE` | Alternate support consistent with `ERR` | `weak_evidence` |
| `NOT_SIGNIFICANT` | Enrichment `EQ` above 0.05 | — |
| `STRAND_BIAS` | `SP` below 0.001, or strand skew above `--max-strand-bias` | `strand_bias` |
| `ORIENTATION_BIAS` | `OP` below 0.001 | `orientation` |
| `BASE_QUAL` / `MAP_QUAL` / `POSITION` | Alternate observations systematically worse than reference on `RSBQ` / `RSMQ` / `RSPOS` | `base_qual`, `map_qual`, `position` |
| `BLACKLIST` | Position in `--custom-blacklist` | — |
| `CIRCULAR_EDGE_UNRESOLVED` | Within `--circular-edge-bases` of a linear reference edge | — |
| `POSSIBLE_NUMT` | Alternate support within `--autosomal-median-depth` | — |

The three rank-sum filters fire on a negative z-score only. A positive skew means
the alternate observations are *better* than the reference ones, which is not
evidence against the allele.

### Input expectations

`paired` does not realign and does not mark duplicates. Give it a BAM or CRAM
that is already deduplicated and, ideally, mate-overlap corrected, and pass
`--deduplication none`. The built-in coordinate deduplication is a fallback for
unmarked input: it is cruder than MarkDuplicates and degenerates to
`(start, strand)` for orphan reads, which over-collapses capture data.

### Scope

`paired` is a mitochondrial heteroplasmy evidence generator, not a general
somatic WES/WGS caller. It deliberately has no local realignment or haplotype
assembly, no indel calling, no panel of normals, no contamination estimate, and
no germline prior. For nuclear somatic calling, use Mutect2, Strelka2, or
DeepSomatic, and read this VCF alongside them.

There is no panel of normals because the arithmetic does not work at
mitochondrial depth. chrM is 16,569 bp, giving 49,707 possible SNV alleles; at
1,000x with a Q30 error rate, 63% of them carry at least one alternate read in
any given normal, rising to nearly all of them by 5,000x. A presence-based PON
would blacklist the whole genome. The workable version is a site-specific
beta-binomial background model, the shearwater/deepSNV approach, which is a
different object altogether and is not implemented.

Other current limitations:

- SNVs only; indels are not emitted.
- Shifted-reference evidence is not yet accepted, so a linear reference leaves
  its artificial breakpoint unresolved. Edge positions stay in the evidence but
  are flagged and excluded from callable territory.
- The bundled NUMT BEDs are nuclear-side masking resources, not mitochondrial
  blacklists. Use `--custom-blacklist` for a chrM-side BED.

### Migrating from paired schema 2.0

| 2.0 | 3.0 |
|---|---|
| `--query` / `--baseline` | `--tumor` / `--normal` |
| `--min-query-depth` / `--min-baseline-depth` | `--min-tumor-depth` / `--min-normal-depth` |
| `--min-query-af` / `--max-baseline-af` | `--min-tumor-af` / `--max-normal-af` |
| `--min-query-baseline-ratio` | removed; it was never applied to any filter |
| `QUERY_*` / `BASELINE_*` columns | `tumor_*` / `normal_*`, lowercase |
| `QUERY` / `BASELINE` VCF samples | `NORMAL`, `TUMOR`, in that order |
| `MBQ` / `MMQ` / `MPOS` in `FORMAT` | `INFO`, `Number=R`, tumour-derived |
| `QBRATIO` | removed; `EQ` already ranks candidates |
| `.mt_candidates.tsv.gz` | removed; every column was already in the VCF |
| `.mt_evidence.tsv.gz` | removed; all-position depth is no longer emitted |
| `.mt_callable.bed.gz` | kept, unchanged |
| `.mt_qc.json` | removed; the same JSON is the `##mgatk2_qc=` header line |
| `.paired.log` | removed; it restated the run summary printed to stdout |

## Contributing and releasing

```bash
make            # format, lint, test, and run the fixtures end to end
python -m build # wheel and source distribution
```

CI repeats Ruff and the tests on Python 3.10–3.12 and builds both distributions.
Tagged releases build and publish the Docker image. Maintainer notes and the
invariants worth knowing before changing anything live in
[AGENTS.md](AGENTS.md).

## Licence

MIT licensed; see [LICENSE](LICENSE). mgatk2 derives from the original
[mgatk](https://github.com/caleblareau/mgatk) by Caleb Lareau (2020).

## Citation

If you use mgatk2, please cite the paper that established the methodology:

> Lareau CA, Ludwig LS, Muus C, et al. Massively parallel single-cell
> mitochondrial DNA genotyping and chromatin profiling. *Nature Biotechnology*
> 39, 451–461 (2021). https://doi.org/10.1038/s41587-020-0645-6

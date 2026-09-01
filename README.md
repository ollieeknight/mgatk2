# mgatk2

mgatk2 is a Python toolkit for mitochondrial genotyping from single-cell and
paired bulk sequencing data. It reimplements the core workflow of
[mgatk](https://github.com/caleblareau/mgatk) without Snakemake or Java and adds
HDF5 output, QC reports, and a tumour/normal SNV evidence workflow.

## Installation

mgatk2 requires Python 3.10 or newer.

```bash
python -m venv .venv
source .venv/bin/activate
pip install git+https://github.com/ollieeknight/mgatk2.git
mgatk2 --help
```

For development:

```bash
git clone https://github.com/ollieeknight/mgatk2.git
cd mgatk2
make setup
make
```

`make` formats, lints, runs the unit tests, and executes every pipeline
end-to-end against the committed fixtures.

## Commands

| Command | Purpose |
|---|---|
| `mgatk2 run` | Single-cell BAM with barcode tags; HDF5 defaults |
| `mgatk2 tenx` | 10x `outs/` input; text defaults compatible with Signac |
| `mgatk2 call` | Directory containing one bulk BAM per sample |
| `mgatk2 paired` | Tumour/normal mitochondrial SNV evidence from BAM or CRAM |
| `mgatk2 hardmask-fasta` | Hard-mask a reference FASTA with bundled NUMT BEDs |

Run `mgatk2 COMMAND --help` for the current options and defaults.

## Single-cell analysis

Use `tenx` for original-mgatk-style text output:

```bash
mgatk2 tenx --input path/to/cellranger-atac/outs --output results
```

Use `run` for filtered HDF5 output and an HTML QC report:

```bash
mgatk2 run \
  --input path/to/possorted_bam.bam \
  --barcodes path/to/barcodes.tsv \
  --output results
```

`--input` also accepts a CellRanger `outs/` directory, for either scATAC
(`possorted_bam.bam`) or 10x Multi/CITE-seq (`per_sample_outs/<sample>/count/
sample_alignments.bam`); the matching barcode file is auto-detected in both
cases. If a Multi `outs/` contains more than one demultiplexed sample, point
`--input` at that sample's `sample_alignments.bam` directly.

If `--barcodes` is omitted, barcodes are extracted from the BAM tag selected by
`--barcode-tag` (default `CB`). The three deduplication modes are:

- `alignment_and_fragment_length`: alignment start, strand, and template length.
- `alignment_start`: alignment start and strand; the `tenx` default.
- `none`: retain all otherwise eligible alignments.

HDF5 output is written below `<output>/output/`:

- `counts.h5`: strand-specific A/C/G/T counts, barcodes, and optional Tn5 cuts.
- `metadata.h5`: coverage, per-cell summaries, inferred reference, and optional
  barcode metadata.

Text output contains compressed base and coverage tables, a depth table, and an
inferred reference allele table. Both formats also write `output.log` and
`qc/{cell_stats.csv,summary.txt}`.

Tn5 cut sites record one insertion per retained read, at the read's outermost
aligned reference base. They are meaningful for ATAC data; disable them with
`--no-tn5` for RNA-seq or WGS.

### How the single-cell scan works

The barcode list is split into contiguous shards. Each shard is handled by one
worker that opens the BAM itself, streams chrM once, and accumulates directly
into dense per-cell count arrays. Reads are never held as Python objects, and a
finished shard travels back as NumPy blocks that land in HDF5 as one contiguous,
chunk-aligned column write.

Two flags control this:

- `--threads` sets the number of concurrent shard workers.
- `--memory` sets the budget those workers share. Shard width is derived from
  it, so peak memory tracks the cell count and the budget, never the number of
  reads in the BAM. A run with more cells than the budget allows simply uses
  more, narrower shards.

Each shard decodes chrM independently, so total decode work scales with the
shard count while wall time does not. That trade buys a single-pass, bounded
memory design with no inter-process read traffic.

The R helpers in [R/mgatk2_functions.R](R/mgatk2_functions.R) load and analyse
HDF5 output. See [R/README.md](R/README.md).

## Paired tumour/normal analysis

`paired` compares a tumour sample against an autologous normal:

```bash
mgatk2 paired \
  --tumor tumour.bam \
  --normal normal.bam \
  --reference GRCh38.fa \
  --output results \
  --sample-name K010
```

The tumour/normal vocabulary is naming, not biology. Any two autologous samples
pair validly — an NKG2C+ population against CD14+ from the same donor is the
same test — but the output uses somatic conventions so it drops straight into
`bcftools`, VEP, and the rest of the somatic tool chain. Column names and VCF
fields use American spelling (`tumor_af`, `normal_dp`) to match that ecosystem.

For upstream UMI-consensus inputs, declare that provenance and do not deduplicate
again:

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

The FASTA defines `REF`. Overlapping mates contribute at most one fragment
observation per position; higher base quality resolves a disagreement and an
equal-quality disagreement is masked. Missing-quality, duplicate, QC-failed,
unmapped, secondary, and supplementary reads are excluded and counted in QC.

Schema 3.0 writes three files and nothing else:

| File | Contents |
|---|---|
| `{sample}.mt_variants.vcf.gz` | Two-sample SNV VCF, the complete variant record |
| `{sample}.mt_variants.vcf.gz.tbi` | Tabix index |
| `{sample}.mt_callable.bed.gz` | Positions reportable in both samples |

The callable BED is the only companion file, because it is the only thing the
VCF cannot stand in for: a VCF carries variant sites, so callable territory is
unrecoverable once the run is over. It is the denominator for any
per-callable-base mutation burden, and the same companion GATK CallableLoci and
Strelka emit alongside a somatic VCF.

### The VCF is the record

Mutect2, Strelka2, and DeepSomatic all emit a single VCF, and so does mgatk2.
Samples are ordered `NORMAL`, `TUMOR`.

- `FORMAT`: `GT:DP:AD:AF:SB:F1R2:F2R1:AFCI`, per sample.
- `INFO`: `DP`, `ERR`, `SEQP`, `EP`, `EQ`, `SP`, `OP`; the tumour-derived
  per-allele medians `MBQ`, `MMQ`, `MPOS` (`Number=R`, reference then
  alternate, matching Mutect2's layout); and the alternate-versus-reference
  rank-sum z-scores and p-values `RSBQ`, `RSMQ`, `RSPOS`.
- `QUAL` is the phred-scaled BH-adjusted enrichment q-value.

Quality medians are tumour-derived only. Nothing filters on the normal-side
values and carrying them doubled every record for no reader.

`AD` carries genuine reference and alternate counts, not `DP - AC`: at a
multiallelic position those differ.

### QC and provenance live in the VCF header

There is no sidecar JSON. The full QC record — parameters, input paths and read
statistics, reference checksum, per-substitution error rates, depth quantiles,
candidate and callable-position counts, mgatk2 version, git commit, and the
command line — is a single `##mgatk2_qc=` header line holding one JSON object.
It survives `bcftools` round-trips:

```bash
bcftools view -h sample.mt_variants.vcf.gz | sed -n 's/^##mgatk2_qc=//p' | jq
```

Loading into R is `VariantAnnotation::readVcf`, or `bcftools query` to a data
frame:

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%DP\t%AD\t%AF]\n' \
  sample.mt_variants.vcf.gz
```

### How a candidate is judged

Two independent questions must both be answered before a variant passes:

1. **Is the tumour enriched over the normal?** One-sided Fisher exact test,
   Benjamini-Hochberg adjusted across every observed allele (`EP`, `EQ`).
2. **Is the tumour alternate above sequencing noise?** A per-substitution error
   rate is learned from the normal sample — excluding sites that carry a
   plausible real allele — and the tumour alternate count is tested against it
   with a binomial (`ERR`, `SEQP`). This matters because a Fisher test alone
   assumes both samples share an error rate, so unequal depth can look
   significant on its own.

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

The three rank-sum filters fire on a negative z-score only: a positive skew
means the alternate observations are *better* than the reference ones, which is
not evidence against the allele.

### Scope

`paired` is a **mitochondrial heteroplasmy evidence generator**, not a general
somatic WES/WGS caller. It deliberately has no local realignment or haplotype
assembly, no indel calling, no panel of normals, no contamination estimation,
and no germline population prior. For nuclear somatic calling use Mutect2,
Strelka2, or DeepSomatic and feed this VCF alongside them.

There is no panel of normals because the arithmetic does not work at
mitochondrial depth. chrM is 16,569 bp, so there are 49,707 possible SNV
alleles; at 1,000x with a Q30 error rate, 63% of them carry at least one
alternate read in any given normal, rising to ~100% by 5,000x. A
presence-based PON would blacklist the entire genome. The workable version is a
site-specific beta-binomial background model rather than a site list — the
shearwater/deepSNV approach — which is a different object from what a PON is,
and is not implemented.

Other current limitations:

- SNVs only; indels are not emitted.
- Shifted-reference evidence is not yet accepted.
- A linear mitochondrial reference leaves its artificial breakpoint unresolved.
  Configured edge positions remain in evidence but are flagged and excluded from
  callable territory.
- Bundled NUMT BEDs are nuclear-side reference-masking resources, not
  mitochondrial blacklists. Use `--custom-blacklist` for a chrM-side BED.

### Input expectations

`paired` does not realign or mark duplicates. Give it a BAM/CRAM that is
already deduplicated and, ideally, mate-overlap corrected, and pass
`--deduplication none`. mgatk2's own coordinate deduplication is a fallback for
unmarked input: it is cruder than MarkDuplicates and degenerates to
`(start, strand)` for orphan reads, which over-collapses capture data.

### Migrating from schema 2.0

| 2.0 | 3.0 |
|---|---|
| `--query` / `--baseline` | `--tumor` / `--normal` |
| `--min-query-depth` / `--min-baseline-depth` | `--min-tumor-depth` / `--min-normal-depth` |
| `--min-query-af` / `--max-baseline-af` | `--min-tumor-af` / `--max-normal-af` |
| `--min-query-baseline-ratio` | removed; it was never applied to any filter |
| `QUERY_*` / `BASELINE_*` columns | `tumor_*` / `normal_*`, lowercase |
| `QUERY` / `BASELINE` VCF samples | `NORMAL`, `TUMOR`, in that order |
| `MBQ`/`MMQ`/`MPOS` in `FORMAT` | `INFO`, `Number=R`, tumour-derived |
| `QBRATIO` | removed; `EQ` already ranks candidates |
| `.mt_candidates.tsv.gz` | removed; every column was already in the VCF |
| `.mt_evidence.tsv.gz` | removed; all-position depth is no longer emitted |
| `.mt_callable.bed.gz` | kept, unchanged |
| `.mt_qc.json` | removed; the same JSON is the `##mgatk2_qc=` header line |
| `.paired.log` | removed; it restated the run summary printed to stdout |

## Releasing

Before tagging:

```bash
make
python -m build
git status --short
```

CI repeats Ruff and the Python tests on 3.10–3.12 and builds the wheel and source
distribution. Tagged releases build and publish the Docker image.

## Licence and citation

mgatk2 is MIT licensed and derives from the original mgatk by Caleb Lareau. If
you use it in research, cite:

> Lareau CA, Ludwig LS, Muus C, et al. Massively parallel single-cell
> mitochondrial DNA genotyping and chromatin profiling. *Nature Biotechnology*
> 39, 451–461 (2021). https://doi.org/10.1038/s41587-020-0645-6

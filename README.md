# mgatk2

mgatk2 is a Python toolkit for mitochondrial genotyping from single-cell and
paired bulk sequencing data. It reimplements the core workflow of
[mgatk](https://github.com/caleblareau/mgatk) without Snakemake or Java and adds
HDF5 output, QC reports, and a query/baseline SNV evidence workflow.

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
| `mgatk2 paired` | Query/baseline mitochondrial SNV evidence from BAM or CRAM |
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

## Paired query/baseline analysis

`paired` compares an experimental query with an autologous baseline without
assuming tumour/germline biology:

```bash
mgatk2 paired \
  --query NKG2C.bam \
  --baseline CD14.bam \
  --reference GRCh38.fa \
  --output results \
  --sample-name K010_NKG2C_v_CD14
```

For upstream UMI-consensus inputs, declare that provenance and do not deduplicate
again:

```bash
mgatk2 paired \
  --query query.consensus.bam \
  --baseline baseline.consensus.bam \
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

Schema 2.0 writes:

| File | Contents |
|---|---|
| `{sample}.mt_evidence.tsv.gz` | A/C/G/T evidence at every mitochondrial position |
| `{sample}.mt_candidates.tsv.gz` | Every observed non-reference SNV, its statistics and filters |
| `{sample}.mt_variants.vcf.gz` and `.tbi` | Indexed two-sample SNV VCF |
| `{sample}.mt_callable.bed.gz` | Positions reportable in both samples |
| `{sample}.mt_qc.json` | Parameters, provenance, schemas, error rates, read/fragment QC |
| `{sample}.paired.log` | Concise execution summary |

The VCF carries everything the candidates table does, using GATK-compatible
field names so it drops straight into `bcftools`, VEP, and downstream union
tooling:

- `FORMAT`: `GT:DP:AD:AF:SB:F1R2:F2R1:MBQ:MMQ:MPOS:AFCI`, per sample. `MBQ`,
  `MMQ`, and `MPOS` are **per-allele** medians of base quality, mapping
  quality, and distance from the nearest read end.
- `INFO`: `DP`, `QBRATIO`, `ERR`, `SEQP`, `EP`, `EQ`, `SP`, `OP`, and the
  alternate-versus-reference rank-sum z-scores and p-values `RSBQ`, `RSMQ`,
  `RSPOS`.
- `QUAL` is the phred-scaled BH-adjusted enrichment q-value.

### How a candidate is judged

Two independent questions must both be answered before a variant passes:

1. **Is the query enriched over the baseline?** One-sided Fisher exact test,
   Benjamini-Hochberg adjusted across every observed allele (`EP`, `EQ`).
2. **Is the query alternate above sequencing noise?** A per-substitution error
   rate is learned from the baseline sample — excluding sites that carry a
   plausible real allele — and the query alternate count is tested against it
   with a binomial (`ERR`, `SEQP`). This matters because a Fisher test alone
   assumes both samples share an error rate, so unequal depth can look
   significant on its own.

Artefact filters then apply: `STRAND_BIAS`, `ORIENTATION_BIAS` (read-pair
F1R2/F2R1 skew, only evaluated when both mates are present), `WEAK_EVIDENCE`,
`BLACKLIST`, `CIRCULAR_EDGE_UNRESOLVED`, and the depth and allele-fraction
thresholds. Supplying `--autosomal-median-depth` additionally enables
`POSSIBLE_NUMT`, which flags alternate support a single-copy NuMT could
account for.

### Scope

`paired` is a **mitochondrial heteroplasmy evidence generator**, not a general
somatic WES/WGS caller. It deliberately has no local realignment or haplotype
assembly, no indel calling, no panel of normals, no contamination estimation,
and no germline population prior. For nuclear somatic calling use Mutect2,
Strelka2, or DeepSomatic and feed this VCF alongside them.

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

The removed `mgatk2 wes` adapter was only a tumour/normal renaming of `paired`.
Replace it with:

```bash
mgatk2 paired --query tumour.cram --baseline normal.cram \
  --reference GRCh38.fa --output results --sample-name pair \
  --deduplication none --input-is-consensus
```

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

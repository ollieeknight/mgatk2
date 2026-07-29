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
make check-all
make integration
```

## Commands

| Command | Purpose |
|---|---|
| `mgatk2 run` | Single-cell BAM with barcode tags; HDF5 defaults |
| `mgatk2 tenx` | 10x `outs/` input; text defaults compatible with Signac |
| `mgatk2 call` | Directory containing one bulk BAM per sample |
| `mgatk2 paired` | Query/baseline mitochondrial SNV evidence from BAM or CRAM |
| `mgatk2 wes` | Deprecated tumour/normal compatibility wrapper for `paired` |
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

Schema 1.0 writes:

| File | Contents |
|---|---|
| `{sample}.mt_evidence.tsv.gz` | A/C/G/T evidence at every mitochondrial position |
| `{sample}.mt_candidates.tsv.gz` | Every observed non-reference SNV and its filters |
| `{sample}.mt_variants.vcf.gz` and `.tbi` | Indexed query/baseline SNV VCF |
| `{sample}.mt_callable.bed.gz` | Positions reportable in both samples |
| `{sample}.mt_qc.json` | Parameters, provenance, schemas, and read/fragment QC |
| `{sample}.paired.log` | Concise execution summary |

The initial numeric thresholds are migration annotations, not validated
biological defaults. `LEGACY_FILTER` records their result; `FILTER` records
technical and statistical flags.

Current paired limitations:

- SNVs only; indels are not emitted.
- Shifted-reference evidence is not yet accepted.
- A linear mitochondrial reference leaves its artificial breakpoint unresolved.
  Configured edge positions remain in evidence but are flagged and excluded from
  callable territory.
- Bundled NUMT BEDs are nuclear-side reference-masking resources, not
  mitochondrial blacklists. Use `--custom-blacklist` for a chrM-side BED.

`mgatk2 wes` remains as a deprecated adapter for one compatibility release. It
maps tumour/normal arguments to the same paired pipeline and additionally writes
`{sample}.mito_somatic.tsv`.

## Releasing

Before tagging:

```bash
make check-all
make integration
git status --short
```

`make check-all` runs Ruff, pytest, CLI smoke tests, and builds both the wheel and
source distribution. CI repeats Python tests on 3.10–3.12. Tagged releases build
and publish the Docker image.

## Licence and citation

mgatk2 is MIT licensed and derives from the original mgatk by Caleb Lareau. If
you use it in research, cite:

> Lareau CA, Ludwig LS, Muus C, et al. Massively parallel single-cell
> mitochondrial DNA genotyping and chromatin profiling. *Nature Biotechnology*
> 39, 451–461 (2021). https://doi.org/10.1038/s41587-020-0645-6

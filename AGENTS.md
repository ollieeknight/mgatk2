# mgatk2 maintainer guide

mgatk2 = Python 3.10+ mitochondrial genotyping toolkit. Public language stay
assay-neutral: paired workflow compare query vs baseline, assume no
tumour/germline biology.

## Layout

```text
src/
  cli/          Click entry point, command modules, and option decorators
  core/         Configuration, exceptions, and single-cell orchestration
  processing/   Sharded chrM scanning, base counting, and fragment grouping
  analysis/     Paired calling, per-allele quality statistics, QC, HTML reports
  file_io/      Sharded single-cell and atomic paired writers
  data/         Bundled nuclear-side NUMT BEDs
  utils/        Input validation and FASTA masking
R/              HDF5 loading and downstream analysis helpers
tests/          Focused pytest suite and committed integration fixtures
```

Generic top-level package names (`cli`, `core`, etc) historic, form current
installed interface. No move without deliberate package migration.

## Required checks

```bash
make setup   # once, creates .venv
make         # format, lint, pytest, end-to-end fixture runs
```

No commit generated outputs; `make` write them under `.test-work/`.

## Single-cell architecture

Single-cell path = sharded streaming scan. Understand before changing anything
in `processing/` or `file_io/writers.py`.

- Barcode list split into **contiguous** shards by `plan_shards`. Shards must
  stay contiguous: writers depend on it to emit chunk-aligned HDF5 column blocks.
- Each shard worker open BAM itself, stream chrM once. Reads never materialised
  as Python objects, never cross process boundary; only finished NumPy blocks
  travel back.
- Shard width derived from `--memory` and `--threads`, so peak memory track cell
  count and budget, never read count. More cells than budget allow = more,
  narrower shards.
- Each shard decode chrM independently. Total decode work scale with shard count,
  wall time not. Deliberate trade for single-pass, bounded-memory design with no
  inter-process read traffic.
- No reintroduce whole-BAM read buffer, per-cell dict payload, batch barrier, or
  `gc.collect()` in shard loop. Each was measured bottleneck.

## Single-cell invariants

- `run` default HDF5 + fragment-length-aware deduplication.
- `tenx` default Signac-compatible text + alignment-start deduplication.
- `call` treat every BAM in its input directory as one bulk sample.
- Deduplication per cell, keyed on alignment start, strand, optionally template
  length.
- CIGAR insertions advance query offset. Lose that = silently shift every base
  after insertion onto wrong reference position.
- Tn5 cut sites record exactly one insertion per retained read, at read's
  outermost aligned reference base. Cut totals must equal retained read count.
  Disable with `--no-tn5`.
- `mean_depth` and `median_depth` average over covered positions only — historic,
  inconsistent with `genome_coverage`. Changing it = QC-visible break.
- HDF5 matrices stored positions × cells; `hdf5r` read transposed cells ×
  positions shape in R.
- Single-cell reference allele inferred from aggregate counts.

## Paired invariants

- FASTA define `REF`; observed major alleles never replace it.
- CRAM decoding require same reference FASTA.
- Mates grouped into fragments, overlapping position count once.
- `paired` no realign, no mark duplicates. Expect already-deduplicated input +
  `--deduplication none`. Built-in coordinate dedup = fallback for unmarked
  input only, cruder than MarkDuplicates, degenerate to `(start, strand)` for
  orphans.
- Evidence contain every FASTA position + raw strand-specific A/C/G/T counts.
- Quality stats stored as per-allele histograms (`analysis/quality_stats.py`),
  never running sums. Pooled ref+alt mean cannot separate real allele from
  artefact — that was the schema 1.0 mistake.
- Candidate pass = two independent tests. `EP`/`EQ` = query enriched over
  baseline. `SEQP` vs `ERR` = query alt above learned substitution error rate.
  Fisher alone assume shared error rate, so depth asymmetry alone look
  significant.
- Error rate learned from baseline, excluding sites with plausible real allele.
  Drop that exclusion = estimate absorb the heteroplasmy it must detect.
- Orientation test only valid when both mates present. Single-end/orphan input
  leave F1R2/F2R1 zero; test must return 1.0, never flag.
- VCF is the primary artefact. GATK-compatible names (`GT:DP:AD:AF:SB:F1R2:
  F2R1:MBQ:MMQ:MPOS:AFCI`), one `##fileformat`, `##reference`, `QUAL` from
  phred-scaled `EQ`. Anything computed and not surfaced in the VCF is a bug.
- TSV and VCF `FILTER` both use `;`.
- Caller SNV-only. Linear-reference edge positions flagged until
  shifted-reference evidence implemented.
- Bundled NUMT BEDs mask nuclear reference regions; only user-supplied chrM BED
  may filter mitochondrial positions. `--autosomal-median-depth` enable
  `POSSIBLE_NUMT`.

Paired schema versions live in `PairedConfig`. Any column or semantic change
need schema-version decision + tests.

Scope: mitochondrial heteroplasmy evidence generator, not general somatic
WES/WGS caller. Deliberately no local realignment, no assembly, no indels, no
panel of normals, no contamination estimate, no germline prior. Do not add
them; feed Mutect2/Strelka2/DeepSomatic instead.

`wes` command removed: was only tumour/normal renaming of `paired` forcing
`--deduplication none --input-is-consensus`. No reinstate assay-specific alias;
add options to `paired` instead.

## Dependencies and release surfaces

Python runtime dependencies declared only in `pyproject.toml`: Click, pysam,
NumPy, tqdm, h5py, Matplotlib, SciPy. Release-relevant files: `README.md`,
`LICENSE`, `pyproject.toml`, `Makefile`, `Dockerfile`, `.dockerignore`, and
workflows under `.github/workflows/`.
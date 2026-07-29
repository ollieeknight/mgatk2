# mgatk2 maintainer guide

mgatk2 is a Python 3.10+ mitochondrial genotyping toolkit. Keep public language
assay-neutral: the paired workflow compares a query with a baseline and does not
assume tumour/germline biology.

## Layout

```text
src/
  cli/          Click entry point, command modules, and option decorators
  core/         Configuration, exceptions, and single-cell orchestration
  processing/   BAM/CRAM readers, fragments, pileups, and cell workers
  analysis/     Paired calling statistics, QC metadata, and HTML reports
  file_io/      Incremental single-cell and atomic paired writers
  data/         Bundled nuclear-side NUMT BEDs
  utils/        Input validation and FASTA masking
R/              HDF5 loading and downstream analysis helpers
tests/          Focused pytest suite and generated integration fixture
```

The generic top-level package names (`cli`, `core`, and so on) are historical and
form the current installed interface. Do not move them without a deliberate
package migration.

## Required checks

```bash
make setup
make check-all
make integration
```

`make check-all` runs formatting checks, Ruff, pytest, CLI smoke tests, and
builds the wheel and sdist. `make integration` generates a tiny ignored 10x-style
fixture and exercises both HDF5 and text pipelines. Do not commit generated
outputs.

## Single-cell invariants

- `run` defaults to HDF5 and fragment-length-aware deduplication.
- `tenx` defaults to Signac-compatible text and alignment-start deduplication.
- `call` treats every BAM in its input directory as one bulk sample.
- BAM reads are streamed from chrM, grouped by barcode, then processed by cell.
- HDF5 matrices are stored as positions × cells; `hdf5r` reads the transposed
  cells × positions shape in R.
- The current single-cell reference allele is inferred from aggregate counts.
- Tn5 cut sites are meaningful for ATAC data and can be disabled with `--no-tn5`.

## Paired invariants

- FASTA defines `REF`; observed major alleles never replace it.
- CRAM decoding requires that same reference FASTA.
- Mates are grouped into fragments and an overlapping position counts once.
- UMI-consensus input requires `--input-is-consensus --deduplication none`.
- Evidence contains every FASTA position and raw strand-specific A/C/G/T counts.
- `LEGACY_FILTER` and technical/statistical `FILTER` are separate.
- The caller is SNV-only. Linear-reference edge positions are flagged until
  shifted-reference evidence is implemented.
- Bundled NUMT BEDs mask nuclear reference regions; only a user-supplied chrM BED
  may filter mitochondrial positions.
- `wes` is a deprecated adapter and must delegate to the paired implementation.

Paired schema versions live in `PairedConfig`. Any column or semantic change
requires a schema-version decision and corresponding tests.

## Dependencies and release surfaces

Python runtime dependencies are declared only in `pyproject.toml`: Click, pysam,
NumPy, tqdm, h5py, Matplotlib, and SciPy. Release-relevant files are `README.md`,
`LICENSE`, `pyproject.toml`, `Makefile`, `Dockerfile`, `.dockerignore`, and the
workflows under `.github/workflows/`.

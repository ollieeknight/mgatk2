"""Writers module for mgatk2 output formats.

Both writers consume whole shards: one contiguous block of barcode columns
arriving as dense arrays. That keeps HDF5 writes chunk-aligned and keeps the
per-position Python work proportional to covered positions only.
"""

import errno
import gzip
import logging
import os
import shutil
import tempfile
import time
from pathlib import Path

import h5py
import numpy as np

from core.config import PipelineConfig

logger = logging.getLogger(__name__)

BASES = ("A", "C", "G", "T")
STRANDS = ("fwd", "rev")


def _cell_stat_rows(result, barcodes) -> list[dict]:
    """QC rows for the cells in a shard that passed the read-count filter."""
    rows = []
    for i in np.flatnonzero(result.kept):
        n_reads = int(result.n_reads[i])
        rows.append(
            {
                "barcode": barcodes[result.offset + i],
                "total_reads": n_reads,
                "total_fragments": n_reads - int(result.n_paired[i]) // 2,
                "mean_depth": float(result.mean_depth[i]),
                "coverage_breadth": float(result.coverage_breadth[i]),
            }
        )
    return rows


class IncrementalHDF5Writer:
    """Writes mgatk results to HDF5, one contiguous shard of columns at a time."""

    # HDF5 can briefly return EAGAIN on networked filesystems.
    MAX_RETRIES = 5
    INITIAL_RETRY_DELAY = 0.1
    MAX_RETRY_DELAY = 5.0
    # Column-chunk width; shard writes span whole chunks, so gzip runs once each.
    CHUNK_CELLS = 128
    CHUNK_CACHE_BYTES = 64 * 1024**2

    def __init__(
        self,
        output_dir: Path,
        config: PipelineConfig,
        barcodes: list[str],
        barcode_metadata=None,
    ):
        self.final_output_dir = output_dir / "output"
        self.final_output_dir.mkdir(exist_ok=True, parents=True)
        self.config = config
        self.barcodes = barcodes
        self.n_barcodes = len(barcodes)
        self.n_positions = config.mito_length
        self.barcode_metadata = barcode_metadata

        staging_parent = Path(os.environ.get("TMPDIR", tempfile.gettempdir()))
        self.staging_dir = Path(tempfile.mkdtemp(prefix="mgatk2_hdf5_", dir=staging_parent))
        self.output_dir = self.staging_dir / "output"
        self.output_dir.mkdir(exist_ok=True, parents=True)

        self.cell_stats: list[dict] = []
        self.base_totals = np.zeros((self.n_positions, 4), dtype=np.int64)
        self.write_error_count = 0
        self.cells_written = 0

        logger.info("Staging HDF5 output in %s", self.staging_dir)
        self._init_hdf5_files()

    def _matrix(self, handle, name):
        handle.create_dataset(
            name,
            shape=(self.n_positions, self.n_barcodes),
            dtype=np.uint16,
            compression="gzip",
            compression_opts=4,
            chunks=(self.n_positions, min(self.CHUNK_CELLS, self.n_barcodes)),
            fillvalue=0,
        )

    def _open_h5(self, name):
        return h5py.File(
            self.output_dir / name,
            "w",
            libver="earliest",
            rdcc_nbytes=self.CHUNK_CACHE_BYTES,
            rdcc_nslots=10007,
        )

    def _init_hdf5_files(self):
        self.counts_file = self._open_h5("counts.h5")
        self.counts_file.attrs["n_cells"] = self.n_barcodes
        self.counts_file.attrs["n_positions"] = self.n_positions
        self.counts_file.attrs["mito_chr"] = self.config.mito_chr
        self.counts_file.create_dataset("barcode", data=np.array(self.barcodes, dtype="S"))

        for base in BASES:
            for strand in STRANDS:
                self._matrix(self.counts_file, f"{base}_{strand}")
        for strand in STRANDS:
            self._matrix(self.counts_file, f"tn5_cuts_{strand}")

        self.metadata_file = self._open_h5("metadata.h5")
        self.metadata_file.attrs["mito_chr"] = self.config.mito_chr
        self.metadata_file.attrs["mito_length"] = self.n_positions
        self._matrix(self.metadata_file, "coverage")
        for name, dtype in (
            ("mean_depth", np.float32),
            ("median_depth", np.float32),
            ("max_depth", np.uint16),
            ("genome_coverage", np.float32),
            ("total_bases", np.float32),
        ):
            self.metadata_file.create_dataset(
                name, shape=(self.n_barcodes,), dtype=dtype, fillvalue=0
            )

    def write_shard(self, result, barcodes):
        """Write one shard's columns; the block is contiguous by construction."""
        lo = result.offset
        hi = lo + result.counts.shape[0]

        def write():
            for b_idx, base in enumerate(BASES):
                for s_idx, strand in enumerate(STRANDS):
                    self.counts_file[f"{base}_{strand}"][:, lo:hi] = np.ascontiguousarray(
                        result.counts[:, :, b_idx, s_idx].T
                    )
            for s_idx, strand in enumerate(STRANDS):
                self.counts_file[f"tn5_cuts_{strand}"][:, lo:hi] = np.ascontiguousarray(
                    result.tn5[:, :, s_idx].T
                )
            self.metadata_file["coverage"][:, lo:hi] = np.ascontiguousarray(result.depth.T)
            self.metadata_file["mean_depth"][lo:hi] = result.mean_depth
            self.metadata_file["median_depth"][lo:hi] = result.median_depth
            self.metadata_file["max_depth"][lo:hi] = result.max_depth
            self.metadata_file["genome_coverage"][lo:hi] = result.coverage_breadth
            self.metadata_file["total_bases"][lo:hi] = result.total_bases

        self._with_retry(write, "shard write")
        self._with_retry(self._flush, "flush")

        self.base_totals += result.base_totals
        self.cell_stats.extend(_cell_stat_rows(result, barcodes))
        self.cells_written += int(np.count_nonzero(result.kept))

    def _flush(self):
        self.counts_file.flush()
        self.metadata_file.flush()

    def _with_retry(self, action, what: str):
        """Retry transient EAGAIN from networked filesystems, then give up."""
        delay = self.INITIAL_RETRY_DELAY
        for attempt in range(self.MAX_RETRIES):
            try:
                action()
                if attempt > 0:
                    logger.info("%s succeeded on attempt %s", what, attempt + 1)
                return
            except OSError as e:
                if e.errno != errno.EAGAIN:
                    raise
                self.write_error_count += 1
                if attempt == self.MAX_RETRIES - 1:
                    logger.error("%s failed after %s attempts", what, self.MAX_RETRIES)
                    raise
                logger.warning(
                    "Temporary %s error (attempt %d/%d), retrying in %.2fs: %s",
                    what,
                    attempt + 1,
                    self.MAX_RETRIES,
                    delay,
                    e,
                )
                time.sleep(delay)
                delay = min(delay * 2, self.MAX_RETRY_DELAY)

    def finalize(self, qc_dir: Path):
        """Write reference alleles and metadata, close, then publish from staging."""
        from .formats import write_cell_stats

        logger.info("Computing reference alleles...")
        best = self.base_totals.argmax(axis=1)
        ref_alleles = np.where(
            self.base_totals.max(axis=1) > 0,
            np.array(BASES, dtype="S1")[best],
            b"N",
        )
        self.metadata_file.create_dataset(
            "reference", data=ref_alleles, compression="gzip", compression_opts=4
        )

        if self.barcode_metadata is not None:
            self._write_barcode_metadata()

        self._with_retry(self._flush, "flush")
        self.counts_file.close()
        self.metadata_file.close()
        self._publish_hdf5_files()

        if self.write_error_count > 0:
            logger.warning(
                "Encountered %d temporary write errors during processing (all recovered)",
                self.write_error_count,
            )

        qc_dir.mkdir(exist_ok=True, parents=True)
        if self.cell_stats:
            write_cell_stats(self.cell_stats, qc_dir / "cell_stats.csv")

    def _write_barcode_metadata(self):
        group = self.metadata_file.create_group("barcode_metadata")
        source_order = {bc: i for i, bc in enumerate(self.barcode_metadata.get("barcode", []))}
        reorder = [source_order[bc] for bc in self.barcodes if bc in source_order]

        for column, values_list in self.barcode_metadata.items():
            try:
                values = np.array([values_list[i] for i in reorder])
                if values.dtype == object or (len(values) and isinstance(values[0], str)):
                    values = values.astype("S")
                group.create_dataset(column, data=values, compression="gzip", compression_opts=4)
            except Exception as e:
                logger.warning("Could not store metadata column '%s': %s", column, e)

    def _publish_hdf5_files(self):
        """Move finalised HDF5 files from local staging to the target output directory."""
        self.final_output_dir.mkdir(exist_ok=True, parents=True)
        for filename in ["counts.h5", "metadata.h5"]:
            dest = self.final_output_dir / filename
            if dest.exists():
                dest.unlink()
            shutil.move(str(self.output_dir / filename), str(dest))
        shutil.rmtree(self.staging_dir, ignore_errors=True)


class IncrementalTextWriter:
    """Writes the original mgatk text format, one shard at a time."""

    def __init__(self, output_dir: Path, config: PipelineConfig, barcodes: list[str]):
        self.output_dir = output_dir / "output"
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.config = config
        self.barcodes = barcodes
        self.cell_stats: list[dict] = []
        self.base_totals = np.zeros((config.mito_length, 4), dtype=np.int64)
        self.cell_depths: dict[str, float] = {}

        self.base_files = {b: open(self.output_dir / f"output.{b}.txt", "w") for b in BASES}
        self.coverage_file = open(self.output_dir / "output.coverage.txt", "w")

    def write_shard(self, result, barcodes):
        for i in range(result.counts.shape[0]):
            barcode = barcodes[result.offset + i]
            self.cell_depths[barcode] = float(result.mean_depth[i])

            depth = result.depth[i]
            covered = np.flatnonzero(depth)
            if covered.size:
                self.coverage_file.write(
                    "".join(f"{p + 1},{barcode},{depth[p]}\n" for p in covered)
                )

            for b_idx, base in enumerate(BASES):
                fwd = result.counts[i, :, b_idx, 0]
                rev = result.counts[i, :, b_idx, 1]
                observed = np.flatnonzero(fwd | rev)
                if observed.size:
                    self.base_files[base].write(
                        "".join(f"{p + 1},{barcode},{fwd[p]},{rev[p]}\n" for p in observed)
                    )

        self.base_totals += result.base_totals
        self.cell_stats.extend(_cell_stat_rows(result, barcodes))

    def finalize(self, qc_dir: Path):
        from .formats import write_cell_stats

        for handle in (*self.base_files.values(), self.coverage_file):
            handle.close()

        logger.info("Compressing output .txt files...")
        for name in [f"output.{b}.txt" for b in BASES] + ["output.coverage.txt"]:
            txt_file = self.output_dir / name
            with (
                open(txt_file, "rb") as f_in,
                gzip.open(self.output_dir / f"{name}.gz", "wb", compresslevel=9) as f_out,
            ):
                shutil.copyfileobj(f_in, f_out)
            txt_file.unlink()

        with open(self.output_dir / "output.depthTable.txt", "w") as f:
            for cell, depth in sorted(self.cell_depths.items()):
                f.write(f"{cell}\t{depth:.2f}\n")

        best = self.base_totals.argmax(axis=1)
        ref_alleles = np.where(self.base_totals.max(axis=1) > 0, np.array(BASES)[best], "N")
        ref_file = self.output_dir / f"{self.config.mito_chr}_refAllele.txt"
        with open(ref_file, "w") as f:
            f.write("pos\tref\n")
            f.write("".join(f"{i + 1}\t{allele}\n" for i, allele in enumerate(ref_alleles)))

        qc_dir.mkdir(exist_ok=True, parents=True)
        if self.cell_stats:
            write_cell_stats(self.cell_stats, qc_dir / "cell_stats.csv")

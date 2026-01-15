"""Writers module for mgatk2 output formats."""

import errno
import gzip
import json
import logging
import shutil
import time
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np

from core.config import PipelineConfig

logger = logging.getLogger(__name__)


class IncrementalHDF5Writer:
    """Writes mgatk results incrementally to HDF5 with batched writes."""

    # Retry configuration for network filesystem operations
    MAX_RETRIES = 5
    INITIAL_RETRY_DELAY = 0.1
    MAX_RETRY_DELAY = 5.0
    BATCH_SIZE = 250

    def __init__(
        self,
        output_dir: Path,
        config: PipelineConfig,
        barcodes: list[str],
        barcode_metadata=None,
    ):
        """Initialise incremental writer with batch buffering."""
        self.output_dir = output_dir / "output"
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.config = config
        self.barcodes = barcodes
        self.barcode_to_idx = {bc: i for i, bc in enumerate(barcodes)}
        self.n_barcodes = len(barcodes)
        self.n_positions = config.mito_length
        self.barcode_metadata = barcode_metadata

        self.cell_stats: list[dict] = []
        self.position_base_counts: dict[int, dict[str, int]] = defaultdict(
            lambda: {"A": 0, "C": 0, "G": 0, "T": 0}
        )
        self.cell_depths: dict[str, int] = {}

        self.counts_file: h5py.File
        self.metadata_file: h5py.File
        self.write_error_count = 0
        self.batch_buffer: list[dict] = []
        self.cells_written = 0

        self._init_hdf5_files()

    def _init_hdf5_files(self):
        """Create HDF5 files with pre-allocated datasets and keep them open."""
        bases = ["A", "C", "G", "T"]
        strands = ["fwd", "rev"]

        file_props = h5py.h5p.create(h5py.h5p.FILE_CREATE)
        file_props.set_userblock(0)

        self.counts_file = h5py.File(
            self.output_dir / "counts.h5",
            "w",
            libver="latest",
        )
        self.counts_file.attrs["n_cells"] = self.n_barcodes
        self.counts_file.attrs["n_positions"] = self.n_positions
        self.counts_file.attrs["mito_chr"] = self.config.mito_chr

        self.counts_file.create_dataset("barcode", data=np.array(self.barcodes, dtype="S"))

        for base in bases:
            for strand in strands:
                self.counts_file.create_dataset(
                    f"{base}_{strand}",
                    shape=(self.n_positions, self.n_barcodes),
                    dtype=np.uint16,
                    compression="gzip",
                    compression_opts=4,
                    chunks=(min(1000, self.n_positions), min(100, self.n_barcodes)),
                    fillvalue=0,
                )

        # Add Tn5 transposition site datasets (Tn5 cut sites)
        for strand in strands:
            self.counts_file.create_dataset(
                f"tn5_cuts_{strand}",
                shape=(self.n_positions, self.n_barcodes),
                dtype=np.uint16,
                compression="gzip",
                compression_opts=4,
                chunks=(min(1000, self.n_positions), min(100, self.n_barcodes)),
                fillvalue=0,
            )

        # Create and keep metadata.h5 open with optimised settings
        self.metadata_file = h5py.File(
            self.output_dir / "metadata.h5",
            "w",
            libver="latest",  # Use latest HDF5 format for better performance
        )
        self.metadata_file.attrs["mito_chr"] = self.config.mito_chr
        self.metadata_file.attrs["mito_length"] = self.n_positions
        self.metadata_file.create_dataset(
            "coverage",
            shape=(self.n_positions, self.n_barcodes),
            dtype=np.uint16,
            compression="gzip",
            compression_opts=4,
            chunks=(min(1000, self.n_positions), min(100, self.n_barcodes)),
            fillvalue=0,
        )
        self.metadata_file.create_dataset(
            "mean_depth", shape=(self.n_barcodes,), dtype=np.float32, fillvalue=0
        )
        self.metadata_file.create_dataset(
            "median_depth", shape=(self.n_barcodes,), dtype=np.float32, fillvalue=0
        )
        self.metadata_file.create_dataset(
            "max_depth", shape=(self.n_barcodes,), dtype=np.uint16, fillvalue=0
        )
        self.metadata_file.create_dataset(
            "genome_coverage", shape=(self.n_barcodes,), dtype=np.float32, fillvalue=0
        )
        self.metadata_file.create_dataset(
            "total_bases", shape=(self.n_barcodes,), dtype=np.float32, fillvalue=0
        )

    def write_cell(self, result: dict):
        """Buffer a single cell's data and write in batches."""
        barcode = result["barcode"]

        if barcode not in self.barcode_to_idx:
            return  # Skip unknown barcodes

        # Add to buffer
        self.batch_buffer.append(result)

        # Store QC if present
        if "qc" in result:
            self.cell_stats.append(result["qc"])

        # Write batch when buffer is full
        if len(self.batch_buffer) >= self.BATCH_SIZE:
            self._write_batch()

    def _write_batch(self):
        """Write accumulated batch of cells to HDF5 files."""
        if not self.batch_buffer:
            return

        bases = ["A", "C", "G", "T"]

        # Pre-allocate arrays for this batch
        batch_size = len(self.batch_buffer)
        batch_indices = []

        # Prepare data structures for batch write
        batch_data = {
            f"{base}_{strand}": np.zeros((self.n_positions, batch_size), dtype=np.uint16)
            for base in bases
            for strand in ["fwd", "rev"]
        }
        batch_coverage = np.zeros((self.n_positions, batch_size), dtype=np.uint16)
        batch_tn5_cuts_fwd = np.zeros((self.n_positions, batch_size), dtype=np.uint16)
        batch_tn5_cuts_rev = np.zeros((self.n_positions, batch_size), dtype=np.uint16)
        batch_mean_depth = np.zeros(batch_size, dtype=np.float32)
        batch_median_depth = np.zeros(batch_size, dtype=np.float32)
        batch_max_depth = np.zeros(batch_size, dtype=np.uint16)
        batch_genome_coverage = np.zeros(batch_size, dtype=np.float32)
        batch_total_bases = np.zeros(batch_size, dtype=np.float32)

        # Process each cell in the batch
        for batch_idx, result in enumerate(self.batch_buffer):
            barcode = result["barcode"]
            pileup = result["pileup"]
            bc_idx = self.barcode_to_idx[barcode]
            batch_indices.append(bc_idx)

            # Calculate statistics
            depths = [counts["depth"] for counts in pileup.values()]
            mean_depth = np.mean(depths) if depths else 0
            median_depth = np.median(depths) if depths else 0
            max_depth = np.max(depths) if depths else 0
            total_bases = np.sum(depths) if depths else 0
            positions_covered = sum(1 for d in depths if d > 0)
            genome_coverage = (
                (positions_covered / self.n_positions * 100) if self.n_positions > 0 else 0
            )
            self.cell_depths[barcode] = mean_depth

            # Fill batch arrays
            for pos, counts in pileup.items():
                pos_idx = pos  # Already 0-based
                pos_1based = pos + 1

                # Coverage
                batch_coverage[pos_idx, batch_idx] = min(counts["depth"], 65535)

                # Tn5 cut sites
                batch_tn5_cuts_fwd[pos_idx, batch_idx] = min(counts.get("tn5_cuts_fwd", 0), 65535)
                batch_tn5_cuts_rev[pos_idx, batch_idx] = min(counts.get("tn5_cuts_rev", 0), 65535)

                # Base counts
                for base in bases:
                    fwd = counts.get(f"{base}_fwd", 0)
                    rev = counts.get(f"{base}_rev", 0)
                    total = counts.get(base, 0)

                    batch_data[f"{base}_fwd"][pos_idx, batch_idx] = min(fwd, 65535)
                    batch_data[f"{base}_rev"][pos_idx, batch_idx] = min(rev, 65535)

                    # Track for reference
                    if total > 0:
                        self.position_base_counts[pos_1based][base] += total

            # Store metadata
            batch_mean_depth[batch_idx] = mean_depth
            batch_median_depth[batch_idx] = median_depth
            batch_max_depth[batch_idx] = min(max_depth, 65535)
            batch_genome_coverage[batch_idx] = genome_coverage
            batch_total_bases[batch_idx] = total_bases

        # Write entire batch to HDF5 (single I/O operation per dataset)
        try:
            for key, data in batch_data.items():
                self._write_batch_with_retry(self.counts_file[key], batch_indices, data)

            self._write_batch_with_retry(
                self.counts_file["tn5_cuts_fwd"], batch_indices, batch_tn5_cuts_fwd
            )
            self._write_batch_with_retry(
                self.counts_file["tn5_cuts_rev"], batch_indices, batch_tn5_cuts_rev
            )
            self._write_batch_with_retry(
                self.metadata_file["coverage"], batch_indices, batch_coverage
            )

            # Write metadata arrays
            for batch_idx, bc_idx in enumerate(batch_indices):
                self.metadata_file["mean_depth"][bc_idx] = batch_mean_depth[batch_idx]
                self.metadata_file["median_depth"][bc_idx] = batch_median_depth[batch_idx]
                self.metadata_file["max_depth"][bc_idx] = batch_max_depth[batch_idx]
                self.metadata_file["genome_coverage"][bc_idx] = batch_genome_coverage[batch_idx]
                self.metadata_file["total_bases"][bc_idx] = batch_total_bases[batch_idx]

            # Flush after each batch
            self._flush_with_retry()

            self.cells_written += len(self.batch_buffer)

        except Exception as e:
            logger.error("Failed to write batch after retries: %s", e)
            raise

        # Clear the buffer
        self.batch_buffer.clear()

    def _write_batch_with_retry(self, dataset, bc_indices: list[int], data: np.ndarray):
        """Write batch of cells to HDF5 dataset with retry logic"""
        delay = self.INITIAL_RETRY_DELAY

        for attempt in range(self.MAX_RETRIES):
            try:
                # Write all columns at once
                for batch_idx, bc_idx in enumerate(bc_indices):
                    dataset[:, bc_idx] = data[:, batch_idx]

                if attempt > 0:
                    logger.info("Batch write succeeded on attempt %s", attempt + 1)
                return
            except OSError as e:
                if e.errno == errno.EAGAIN or e.errno == 11:  # Resource temporarily unavailable
                    self.write_error_count += 1
                    if attempt < self.MAX_RETRIES - 1:
                        logger.warning(
                            "Temporary batch write error (attempt %d/%d), retrying in %.2fs: %s",
                            attempt + 1,
                            self.MAX_RETRIES,
                            delay,
                            e,
                        )
                        time.sleep(delay)
                        delay = min(delay * 2, self.MAX_RETRY_DELAY)  # Exponential backoff
                    else:
                        logger.error(f"Batch write failed after {self.MAX_RETRIES} attempts")
                        raise
                else:
                    # Different error, don't retry
                    raise

    def _flush_with_retry(self):
        """Flush HDF5 files with retry logic."""
        delay = self.INITIAL_RETRY_DELAY

        for attempt in range(self.MAX_RETRIES):
            try:
                self.counts_file.flush()
                self.metadata_file.flush()
                return
            except OSError as e:
                if e.errno == errno.EAGAIN or e.errno == 11:
                    if attempt < self.MAX_RETRIES - 1:
                        logger.warning(
                            "Flush error (attempt %d/%d), retrying in %.2fs",
                            attempt + 1,
                            self.MAX_RETRIES,
                            delay,
                        )
                        time.sleep(delay)
                        delay = min(delay * 2, self.MAX_RETRY_DELAY)
                    else:
                        logger.error("Flush failed after %s attempts", self.MAX_RETRIES)
                        raise
                else:
                    raise

    def finalize(self, qc_dir: Path):
        """
        Finalize HDF5 files.

        Flush any remaining batch, write reference alleles, close files, and write QC.

        Args:
            qc_dir: QC output directory
        """
        from .formats import write_cell_stats

        # Write any remaining buffered cells
        if self.batch_buffer:
            self._write_batch()

        # Determine reference alleles more efficiently
        logger.info("Computing reference alleles...")
        bases = ["A", "C", "G", "T"]
        ref_alleles = ["N"] * self.n_positions

        for pos, counts in self.position_base_counts.items():
            if 1 <= pos <= self.n_positions:
                max_base = max(bases, key=lambda b: counts[b])
                if counts[max_base] > 0:
                    ref_alleles[pos - 1] = max_base

        # Write reference to metadata.h5
        ref_array = np.array(ref_alleles, dtype="S1")
        self.metadata_file.create_dataset(
            "reference", data=ref_array, compression="gzip", compression_opts=4
        )

        # Write barcode metadata if available (from singlecell.csv)
        if self.barcode_metadata is not None:
            # Create a group for barcode metadata
            if "barcode_metadata" in self.metadata_file:
                del self.metadata_file["barcode_metadata"]
            metadata_group = self.metadata_file.create_group("barcode_metadata")

            # metadata is a dict with column names as keys and lists as values
            # Create a barcode-to-index mapping for reordering
            barcode_list = self.barcode_metadata.get("barcode", [])
            barcode_to_idx = {bc: i for i, bc in enumerate(barcode_list)}

            # Create index mapping to reorder metadata to match self.barcodes order
            reorder_indices = [barcode_to_idx[bc] for bc in self.barcodes if bc in barcode_to_idx]

            # Store each column from the metadata dict
            for col, values_list in self.barcode_metadata.items():
                try:
                    # Reorder to match barcode order
                    reordered_values = [values_list[i] for i in reorder_indices]
                    values = np.array(reordered_values)

                    # Handle different data types
                    if values.dtype == object or (len(values) > 0 and isinstance(values[0], str)):
                        # Convert strings to bytes for HDF5
                        values = np.array(reordered_values, dtype="S")

                    metadata_group.create_dataset(
                        col, data=values, compression="gzip", compression_opts=4
                    )
                except Exception as e:
                    logger.warning("Could not store metadata column '%s': %s", col, e)

        # Final flush before closing
        self._flush_with_retry()

        # Close HDF5 files
        self.counts_file.close()
        self.metadata_file.close()

        if self.write_error_count > 0:
            logger.warning(
                "Encountered %d temporary write errors during processing (all recovered)",
                self.write_error_count,
            )

        # Write QC files
        qc_dir.mkdir(exist_ok=True, parents=True)
        if self.cell_stats:
            write_cell_stats(self.cell_stats, qc_dir / "cell_stats.csv")


class IncrementalTextWriter:
    """Writes original mgatk format incrementally."""

    def __init__(self, output_dir: Path, config: PipelineConfig, barcodes: list[str]):
        self.output_dir = output_dir / "output"
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.config = config
        self.barcodes = set(barcodes)
        self.cell_stats: list[dict] = []
        self.position_base_counts: dict[int, dict[str, int]] = defaultdict(
            lambda: {"A": 0, "C": 0, "G": 0, "T": 0}
        )
        self.cell_depths: dict[str, int] = {}

        bases = ["A", "C", "G", "T"]
        self.base_files = {}
        for base in bases:
            self.base_files[base] = open(self.output_dir / f"{base}.txt", "w")

        self.coverage_file = open(self.output_dir / "coverage.txt", "w")

    def write_cell(self, result: dict):
        barcode = result["barcode"]
        pileup = result["pileup"]

        if "qc" in result:
            self.cell_stats.append(result["qc"])

        depths = [counts["depth"] for counts in pileup.values()]
        self.cell_depths[barcode] = np.mean(depths) if depths else 0

        pos_data: dict = defaultdict(lambda: defaultdict(lambda: [0, 0]))

        for pos, counts in pileup.items():
            pos_1based = pos + 1
            total_cov = counts["depth"]

            if total_cov > 0:
                self.coverage_file.write(f"{pos_1based},{barcode},{total_cov}\n")

            for base in ["A", "C", "G", "T"]:
                fwd = counts.get(f"{base}_fwd", 0)
                rev = counts.get(f"{base}_rev", 0)
                total = counts.get(base, 0)

                if fwd > 0 or rev > 0:
                    pos_data[base][pos_1based] = [fwd, rev]

                if total > 0:
                    self.position_base_counts[pos_1based][base] += total

        for base, positions in pos_data.items():
            for pos, (fwd, rev) in positions.items():
                self.base_files[base].write(f"{pos},{barcode},{fwd},{rev}\n")

    def finalize(self, qc_dir: Path):
        from .formats import write_cell_stats

        for f in self.base_files.values():
            f.close()
        self.coverage_file.close()

        logger.info("Compressing output files...")
        for base in ["A", "C", "G", "T"]:
            txt_file = self.output_dir / f"{base}.txt"
            gz_file = self.output_dir / f"output.{base}.txt.gz"
            with open(txt_file, "rb") as f_in:
                with gzip.open(gz_file, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            txt_file.unlink()

        txt_file = self.output_dir / "coverage.txt"
        gz_file = self.output_dir / "output.coverage.txt.gz"
        with open(txt_file, "rb") as f_in:
            with gzip.open(gz_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        txt_file.unlink()

        depth_file = self.output_dir / "output.depthTable.txt"
        with open(depth_file, "w") as f:
            for cell, depth in sorted(self.cell_depths.items()):
                f.write(f"{cell}\t{depth:.2f}\n")

        ref_alleles = {}
        for pos in range(1, self.config.mito_length + 1):
            if pos in self.position_base_counts:
                counts = self.position_base_counts[pos]
                max_base = max(["A", "C", "G", "T"], key=lambda b: counts[b])
                ref_alleles[pos] = max_base if counts[max_base] > 0 else "N"
            else:
                ref_alleles[pos] = "N"

        ref_file = self.output_dir / f"{self.config.mito_chr}_refAllele.txt"
        with open(ref_file, "w") as f:
            f.write("pos\tref\n")
            for pos in range(1, self.config.mito_length + 1):
                f.write(f"{pos}\t{ref_alleles[pos]}\n")

        qc_dir.mkdir(exist_ok=True, parents=True)
        if self.cell_stats:
            write_cell_stats(self.cell_stats, qc_dir / "cell_stats.csv")


def write_parameters_json(parameters: dict, output_path: Path):
    """Write analysis parameters to JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump(parameters, f, indent=2, default=str)

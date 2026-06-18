"""WES somatic mitochondrial variant calling — pileup, comparison, and output."""

from __future__ import annotations

import csv
import logging
import time
from pathlib import Path

import numpy as np
import pysam
from tqdm import tqdm

from core.config import SimpleRead, WesConfig
from core.exceptions import BAMReadError, InvalidInputError, NoChrMReadsError

logger = logging.getLogger(__name__)

# TSV output column order
TSV_COLUMNS = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "NORMAL_DEPTH",
    "TUMOUR_DEPTH",
    "NORMAL_REF_COUNT",
    "NORMAL_ALT_COUNT",
    "NORMAL_FWD",
    "NORMAL_REV",
    "TUMOUR_REF_COUNT",
    "TUMOUR_ALT_COUNT",
    "TUMOUR_FWD",
    "TUMOUR_REV",
    "NORMAL_VAF",
    "TUMOUR_VAF",
    "TN_RATIO",
    "STRAND_BIAS",
    "FILTER",
]


class CRAMReader:
    """Reads all chrM reads from a CRAM or BAM file (no barcode filtering)."""

    def __init__(self, cram_path: str, config: WesConfig, reference: str):
        self.cram_path = Path(cram_path)
        self.config = config
        self.reference = reference

        if not self.cram_path.exists():
            raise InvalidInputError(f"Input file not found: {cram_path}")

        ref_path = Path(reference)
        if not ref_path.exists():
            raise InvalidInputError(f"Reference FASTA not found: {reference}")

        # Check index (.crai for CRAM, .bai for BAM)
        is_cram = str(self.cram_path).endswith(".cram")
        index_path = Path(str(self.cram_path) + (".crai" if is_cram else ".bai"))
        if not index_path.exists():
            raise InvalidInputError(
                f"Index not found: {index_path}\nRun: samtools index {cram_path}"
            )

        self._is_cram = is_cram
        self._validate_and_detect_mito()

    def _open(self):
        """Open the alignment file with correct mode."""
        if self._is_cram:
            return pysam.AlignmentFile(str(self.cram_path), "rc", reference_filename=self.reference)
        return pysam.AlignmentFile(str(self.cram_path), "rb")

    def _validate_and_detect_mito(self):
        try:
            aln = self._open()
        except Exception as e:
            raise BAMReadError(str(self.cram_path), f"Cannot open: {e}") from e

        available = list(aln.references)
        for name in [self.config.mito_chr, "chrM", "MT", "M", "chrMT"]:
            if name in available:
                if self.config.mito_chr != name:
                    logger.info("Using mitochondrial chromosome: %s", name)
                    self.config.mito_chr = name
                break
        else:
            aln.close()
            raise NoChrMReadsError(str(self.cram_path), available)

        aln.close()

    def collect_reads(self, label: str = "reads") -> tuple[list[SimpleRead], dict]:
        """Collect all chrM reads. Returns (reads, stats_dict)."""
        reads: list[SimpleRead] = []
        total = 0
        kept = 0
        skipped_mapq = 0

        try:
            with self._open() as aln:
                for read in tqdm(
                    aln.fetch(self.config.mito_chr),
                    desc=f"Reading {label}",
                    unit=" read",
                    bar_format="{desc}: {n:,} reads [{elapsed}, {rate_fmt}]",
                ):
                    total += 1

                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    if read.query_sequence is None:
                        continue
                    if read.mapping_quality < self.config.min_mapq:
                        skipped_mapq += 1
                        continue

                    qualities = (
                        np.array(read.query_qualities, dtype=np.int8)
                        if read.query_qualities is not None
                        else np.full(len(read.query_sequence), 60, dtype=np.int8)
                    )

                    reads.append(
                        SimpleRead(
                            reference_start=read.reference_start,
                            is_reverse=read.is_reverse,
                            mapping_quality=read.mapping_quality,
                            query_sequence=read.query_sequence.encode("ascii"),
                            query_qualities=qualities,
                            cigar=read.cigartuples if read.cigartuples else [],
                            is_proper_pair=read.is_proper_pair,
                            is_paired=read.is_paired,
                            template_length=read.template_length or 0,
                        )
                    )
                    kept += 1

        except Exception as e:
            raise BAMReadError(str(self.cram_path), f"Read error: {e}") from e

        logger.info(
            "%s: %d / %d reads kept on %s (%d low-MAPQ filtered)",
            label,
            kept,
            total,
            self.config.mito_chr,
            skipped_mapq,
        )
        return reads, {"total": total, "kept": kept, "skipped_mapq": skipped_mapq}


def infer_reference_alleles(
    tumour_pileup: dict[int, dict],
    normal_pileup: dict[int, dict],
    mito_length: int,
) -> list[str]:
    """Infer reference allele per position as majority base across both pileups.

    Returns a list of length mito_length with 0-based indexing ('N' if no coverage).
    """
    ref: list[str] = ["N"] * mito_length
    for pos in range(mito_length):
        counts: dict[str, int] = {"A": 0, "C": 0, "G": 0, "T": 0}
        for pileup in (tumour_pileup, normal_pileup):
            p = pileup.get(pos, {})
            for base in "ACGT":
                counts[base] += p.get(base, 0)
        best = max(counts, key=counts.__getitem__)
        if counts[best] > 0:
            ref[pos] = best
    return ref


def call_somatic_variants(
    tumour_pileup: dict[int, dict],
    normal_pileup: dict[int, dict],
    ref_alleles: list[str],
    config: WesConfig,
    blacklist: set[int],
) -> list[dict]:
    """Compare tumour vs normal pileups and return all candidate variant rows.

    Returns ALL candidate sites (including filtered ones) with FILTER column set.
    PASS sites pass all thresholds. Filtered sites have pipe-delimited reasons.
    """
    results = []

    for pos_0 in range(config.mito_length):
        pos_1 = pos_0 + 1  # 1-based for output
        ref = ref_alleles[pos_0]
        if ref == "N":
            continue

        t = tumour_pileup.get(pos_0, {})
        n = normal_pileup.get(pos_0, {})
        t_depth = t.get("depth", 0)
        n_depth = n.get("depth", 0)

        for alt in "ACGT":
            if alt == ref:
                continue

            t_fwd = t.get(f"{alt}_fwd", 0)
            t_rev = t.get(f"{alt}_rev", 0)
            t_alt = t_fwd + t_rev
            n_fwd = n.get(f"{alt}_fwd", 0)
            n_rev = n.get(f"{alt}_rev", 0)
            n_alt = n_fwd + n_rev

            # Skip if zero evidence in both samples (saves rows in output)
            if t_alt == 0 and n_alt == 0:
                continue

            t_ref = t_depth - t_alt
            n_ref = n_depth - n_alt

            t_vaf = t_alt / t_depth if t_depth > 0 else 0.0
            n_vaf = n_alt / n_depth if n_depth > 0 else 0.0
            tn_ratio = t_vaf / n_vaf if n_vaf > 0 else float("inf")

            # Tumour strand bias: |fwd-rev| / (fwd+rev)
            t_bias = abs(t_fwd - t_rev) / (t_fwd + t_rev) if (t_fwd + t_rev) > 0 else 0.0

            # Collect filter failures
            filters = []
            if t_depth < config.min_tumour_depth:
                filters.append("LOW_TUMOUR_DEPTH")
            if n_depth < config.min_normal_depth:
                filters.append("LOW_NORMAL_DEPTH")
            if t_alt < config.min_tumour_alt_reads:
                filters.append("LOW_ALT_READS")
            if t_vaf < config.min_tumour_af:
                filters.append("LOW_VAF")
            if n_vaf > config.max_normal_af:
                filters.append("GERMLINE")
            if tn_ratio != float("inf") and tn_ratio < config.min_tn_ratio:
                filters.append("LOW_TN_RATIO")
            if t_bias > config.max_strand_bias:
                filters.append("STRAND_BIAS")
            if pos_1 in blacklist:
                filters.append("BLACKLIST")

            results.append(
                {
                    "CHROM": config.mito_chr,
                    "POS": pos_1,
                    "REF": ref,
                    "ALT": alt,
                    "NORMAL_DEPTH": n_depth,
                    "TUMOUR_DEPTH": t_depth,
                    "NORMAL_REF_COUNT": n_ref,
                    "NORMAL_ALT_COUNT": n_alt,
                    "NORMAL_FWD": n_fwd,
                    "NORMAL_REV": n_rev,
                    "TUMOUR_REF_COUNT": t_ref,
                    "TUMOUR_ALT_COUNT": t_alt,
                    "TUMOUR_FWD": t_fwd,
                    "TUMOUR_REV": t_rev,
                    "NORMAL_VAF": round(n_vaf, 6),
                    "TUMOUR_VAF": round(t_vaf, 6),
                    "TN_RATIO": round(tn_ratio, 3) if tn_ratio != float("inf") else "Inf",
                    "STRAND_BIAS": round(t_bias, 4),
                    "FILTER": "PASS" if not filters else "|".join(filters),
                }
            )

    return results


def write_variants_tsv(variants: list[dict], output_path: Path) -> int:
    """Write variant list to TSV. Returns number of PASS variants written."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=TSV_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(variants)
    return sum(1 for v in variants if v["FILTER"] == "PASS")


def run_wes_pipeline(
    tumour_cram: str,
    normal_cram: str,
    reference: str,
    output_dir: str,
    config: WesConfig,
    sample_name: str = "sample",
) -> dict:
    """End-to-end WES mito somatic calling pipeline.

    Returns dict with summary stats and output file path.
    """
    from data.blacklists import load_blacklist_positions
    from processing.pileup import PileupGenerator

    t0 = time.time()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load blacklist (empty by default — nuclear NuMT BEDs have no chrM rows)
    logger.info("Loading blacklist (build=%s)...", config.blacklist_build)
    blacklist = load_blacklist_positions(
        build=config.blacklist_build,
        custom_bed=config.custom_blacklist,
        mito_chr=config.mito_chr,
    )
    logger.info("  %d blacklisted chrM positions", len(blacklist))

    pipeline_config = config.to_pipeline_config()
    pileup_gen = PileupGenerator(pipeline_config)

    # Read and pileup tumour
    tumour_reader = CRAMReader(tumour_cram, config, reference)
    tumour_reads, tumour_stats = tumour_reader.collect_reads(label="tumour")
    logger.info("Generating tumour pileup...")
    tumour_pileup = pileup_gen.generate_pileup(tumour_reads)
    del tumour_reads

    # Read and pileup normal
    normal_reader = CRAMReader(normal_cram, config, reference)
    normal_reads, normal_stats = normal_reader.collect_reads(label="normal")
    logger.info("Generating normal pileup...")
    normal_pileup = pileup_gen.generate_pileup(normal_reads)
    del normal_reads

    # Infer reference alleles from majority base in both pileups
    ref_alleles = infer_reference_alleles(tumour_pileup, normal_pileup, config.mito_length)
    n_ref_covered = sum(1 for r in ref_alleles if r != "N")
    logger.info(
        "Reference alleles inferred for %d / %d positions", n_ref_covered, config.mito_length
    )

    # Call variants
    logger.info("Calling somatic variants...")
    variants = call_somatic_variants(tumour_pileup, normal_pileup, ref_alleles, config, blacklist)
    total = len(variants)
    n_pass = sum(1 for v in variants if v["FILTER"] == "PASS")
    logger.info("  %d candidate sites, %d PASS", total, n_pass)

    # Write output
    tsv_path = output_path / f"{sample_name}.mito_somatic.tsv"
    write_variants_tsv(variants, tsv_path)
    logger.info("Written: %s", tsv_path)

    elapsed = time.time() - t0
    logger.info("Done in %.1fs", elapsed)

    return {
        "total_candidates": total,
        "pass_variants": n_pass,
        "tumour_reads_kept": tumour_stats["kept"],
        "normal_reads_kept": normal_stats["kept"],
        "positions_with_ref": n_ref_covered,
        "output_tsv": str(tsv_path),
        "elapsed_s": round(elapsed, 1),
    }

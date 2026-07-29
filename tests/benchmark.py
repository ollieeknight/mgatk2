#!/usr/bin/env python3
"""Side-by-side benchmark: classic vs fast pileup on the test dataset.

Usage:
    cd tests && python benchmark.py [--format txt|hdf5] [--runs N]

Runs both pileup modes, times them, measures peak RSS, then diffs output
to confirm correctness. Prints a summary table.
"""

import argparse
import gzip
import platform
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
REPO_ROOT = SCRIPT_DIR.parent
# Prefer the project environment, then fall back to PATH.
_venv_bin = REPO_ROOT / ".venv" / "bin" / "mgatk2"
MGATK2 = _venv_bin if _venv_bin.exists() else Path(shutil.which("mgatk2") or "mgatk2")
OUTS = REPO_ROOT / ".test-work" / "outs"
TXT_FILES = [
    "output.A.txt.gz",
    "output.C.txt.gz",
    "output.G.txt.gz",
    "output.T.txt.gz",
    "output.coverage.txt.gz",
]


def run_mode(mode: str, out_dir: Path, fmt: str, extra_args: list[str]) -> tuple[float, int]:
    """Run mgatk2 run with given pileup mode. Returns (wall_seconds, peak_rss_mb)."""
    if out_dir.exists():
        shutil.rmtree(out_dir)

    cmd = [
        str(MGATK2),
        "run",
        "-i",
        str(OUTS),
        "-o",
        str(out_dir),
        "-f",
        fmt,
        "--pileup-mode",
        mode,
        *extra_args,
    ]

    t0 = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - t0

    if result.returncode != 0:
        print(f"\n[{mode}] FAILED (exit {result.returncode})")
        print(result.stderr[-2000:])
        sys.exit(1)

    # peak RSS via /proc/self/status is not available for subprocess directly;
    # use resource module on the last rusage — not perfect but indicative on macOS.
    # On Linux we can read the child's /proc/<pid>/status peak_rss after the fact.
    rss_mb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if platform.system() == "Linux":
        rss_mb = rss_mb // 1024  # Linux: kB → MB
    else:
        rss_mb = rss_mb // (1024 * 1024)  # macOS: bytes → MB

    return wall, rss_mb


def parse_txt_gz(path: Path) -> dict[tuple, int]:
    """Parse pos,barcode,fwd,rev into {(pos,barcode,strand): count}."""
    data = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.strip().split(",")
            if len(parts) != 4:
                continue
            pos, bc, fwd, rev = parts
            data[(pos, bc, "fwd")] = int(fwd)
            data[(pos, bc, "rev")] = int(rev)
    return data


def diff_txt_outputs(dir_a: Path, dir_b: Path) -> tuple[bool, int, int]:
    """Compare txt output files between two run dirs.
    Returns (match, n_mismatches, max_diff).
    """
    out_a = dir_a / "output"
    out_b = dir_b / "output"
    total_mismatches = 0
    max_diff = 0

    for fname in TXT_FILES:
        fa = out_a / fname
        fb = out_b / fname
        if not fa.exists() or not fb.exists():
            return False, -1, -1

        da = parse_txt_gz(fa)
        db = parse_txt_gz(fb)

        all_keys = set(da) | set(db)
        for k in all_keys:
            va = da.get(k, 0)
            vb = db.get(k, 0)
            diff = abs(va - vb)
            if diff > 0:
                total_mismatches += 1
                max_diff = max(max_diff, diff)

    return total_mismatches == 0, total_mismatches, max_diff


def count_cells(out_dir: Path) -> int:
    qc = out_dir / "qc" / "cell_stats.csv"
    if not qc.exists():
        return -1
    with open(qc) as f:
        return sum(1 for line in f) - 1  # subtract header


def print_table(rows: list[dict]):
    labels = [
        ("Pileup mode", "mode"),
        ("Wall time (s)", "wall"),
        ("Peak RSS (MB)", "rss"),
        ("Cells processed", "cells"),
        ("Output match", "match"),
        ("Max count diff", "max_diff"),
        ("Mismatches", "mismatches"),
    ]
    col_w = max(len(label) for label, _ in labels) + 2
    val_w = 14

    sep = "+" + "-" * (col_w + 2) + "+" + (("+" + "-" * (val_w + 2)) * len(rows)) + "+"
    print("\n" + sep)
    header = f"| {'Metric':<{col_w}} |"
    for r in rows:
        header += f" {r['mode']:^{val_w}} |"
    print(header)
    print(sep)

    for label, key in labels:
        row = f"| {label:<{col_w}} |"
        for r in rows:
            v = r.get(key, "—")
            if isinstance(v, float):
                v = f"{v:.2f}"
            row += f" {str(v):^{val_w}} |"
        print(row)

    print(sep)

    # speedup
    if len(rows) == 2 and rows[0]["wall"] and rows[1]["wall"]:
        speedup = rows[0]["wall"] / rows[1]["wall"]
        print(f"\n  Speedup (classic → fast): {speedup:.2f}×")


def main():
    parser = argparse.ArgumentParser(description="mgatk2 pileup benchmark")
    parser.add_argument("--format", default="txt", choices=["txt", "hdf5"])
    parser.add_argument("--runs", type=int, default=1, help="Repetitions per mode (best of N)")
    parser.add_argument("--extra", nargs="*", default=[], help="Extra mgatk2 args")
    args = parser.parse_args()

    if not MGATK2.exists() and not shutil.which("mgatk2"):
        sys.exit("mgatk2 not found — activate the environment or run `make setup`")
    if not OUTS.exists():
        sys.exit(f"Test data not found at {OUTS}")

    print("\nmgatk2 pileup benchmark")
    print(
        f"  Platform: {platform.system()} {platform.machine()} | Python {platform.python_version()}"
    )
    print(f"  Format:   {args.format}")
    print(f"  Runs:     {args.runs}")
    print(f"  Extra:    {args.extra or 'none'}")

    modes = ["classic", "fast"]
    results = []
    out_dirs = {}

    for mode in modes:
        best_wall = float("inf")
        best_rss = 0
        out_dir = SCRIPT_DIR / f"benchmark_{mode}_output"
        out_dirs[mode] = out_dir

        for run_i in range(args.runs):
            print(f"\n  [{mode}] run {run_i + 1}/{args.runs} ...", flush=True)
            wall, rss = run_mode(mode, out_dir, args.format, args.extra)
            print(f"    wall={wall:.2f}s  rss={rss}MB")
            if wall < best_wall:
                best_wall = wall
                best_rss = rss

        cells = count_cells(out_dir)
        results.append(
            {
                "mode": mode,
                "wall": best_wall,
                "rss": best_rss,
                "cells": cells,
                "match": "—",
                "max_diff": "—",
                "mismatches": "—",
            }
        )

    # diff outputs (txt mode only; hdf5 diff is more complex)
    if args.format == "txt":
        print("\n  Diffing outputs ...", flush=True)
        match, n_mis, max_diff = diff_txt_outputs(out_dirs["classic"], out_dirs["fast"])
        results[1]["match"] = "YES ✓" if match else "DIFFER ✗"
        results[1]["max_diff"] = max_diff if not match else 0
        results[1]["mismatches"] = n_mis if not match else 0
        results[0]["match"] = "(reference)"
    else:
        print("  (HDF5 diff skipped — compare manually with h5diff)")

    print_table(results)

    # clean up benchmark outputs
    for d in out_dirs.values():
        if d.exists():
            shutil.rmtree(d)

    if args.format == "txt" and results[1]["match"] not in ("YES ✓", "—"):
        sys.exit(1)


if __name__ == "__main__":
    main()

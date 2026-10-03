"""HTML QC report generation for single-cell runs."""

import base64
import html
import json
import logging
import os
import sys
import tempfile
from datetime import datetime
from io import BytesIO
from pathlib import Path

import h5py
import numpy as np

_cache_root = Path(tempfile.gettempdir()) / "mgatk2-cache"
os.environ.setdefault("MPLCONFIGDIR", str(_cache_root / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_cache_root))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter  # noqa: E402

logging.getLogger("matplotlib").setLevel(logging.WARNING)

logger = logging.getLogger(__name__)

# Every plot is inlined as base64, so resolution sets the report's file size.
# 150 dpi is indistinguishable on screen and roughly a quarter of 300 dpi.
PLOT_DPI = 150

_DEPTH_TICKS = [1, 5, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500, 1000, 5000, 10000]


def plot_to_base64(fig):
    """Render a figure as an inline PNG data URI and close it."""
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def _tidy(ax, xlabel, ylabel):
    """Label the axes and apply the shared styling."""
    ax.set_xlabel(xlabel, fontsize=10, color="black")
    ax.set_ylabel(ylabel, fontsize=10, color="black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _empty_plot(message, figsize=(6, 5)):
    fig, ax = plt.subplots(figsize=figsize)
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=14, color="gray")
    ax.axis("off")
    return plot_to_base64(fig)


def _track_plot(ylabel, draw, length):
    """A per-position track along chrM; `draw(ax, positions)` adds the data."""
    positions = np.arange(1, length + 1)
    fig, ax = plt.subplots(figsize=(10, 3.5))
    draw(ax, positions)
    ax.set_xlim(0, length)
    _tidy(ax, "chrM (bp)", ylabel)
    return plot_to_base64(fig)


def _sum_over_cells(dataset) -> np.ndarray:
    """Total a positions x cells matrix per position, one column block at a time.

    Reading the whole matrix costs cells x positions x 2 bytes, which is
    hundreds of megabytes on a real run and was paid once per plot.
    """
    n_positions, n_cells = dataset.shape
    chunk_cells = (dataset.chunks or (n_positions, 1024))[1]
    # Read whole chunks, but several at a time: the column chunks are 128 cells
    # wide, and one read per chunk is all Python overhead for no benefit.
    block = chunk_cells * max(1, 1024 // chunk_cells)
    totals = np.zeros(n_positions, dtype=np.int64)
    for start in range(0, n_cells, block):
        totals += dataset[:, start : start + block].sum(axis=1, dtype=np.int64)
    return totals


def load_report_data(counts_file: Path, metadata_file: Path) -> dict:
    """Read every array the report needs, opening each HDF5 file once."""
    data: dict = {}

    with h5py.File(counts_file, "r") as f:
        data["tn5_fwd"] = _sum_over_cells(f["tn5_cuts_fwd"])
        data["tn5_rev"] = _sum_over_cells(f["tn5_cuts_rev"])

    with h5py.File(metadata_file, "r") as f:
        data["coverage_total"] = _sum_over_cells(f["coverage"])
        data["mean_depth"] = f["mean_depth"][:]
        data["genome_coverage"] = f["genome_coverage"][:]
        data["total_bases"] = f["total_bases"][:]
        data["reference"] = [base.decode() for base in f["reference"][:]]
        data["total_fragments"] = (
            f["barcode_metadata"]["total"][:] if "barcode_metadata" in f else None
        )

    data["n_cells"] = len(data["mean_depth"])
    data["coverage_mean"] = data["coverage_total"] / max(data["n_cells"], 1)
    return data


def create_coverage_plot(mean_coverage):
    """Mean depth along chrM."""

    def draw(ax, positions):
        ax.plot(positions, mean_coverage, linewidth=0.8, color="#2E86AB")

    return _track_plot("Mean depth", draw, len(mean_coverage))


def create_transposition_frequency_plot(tn5_fwd, tn5_rev):
    """Mirrored per-strand Tn5 cut sites along chrM."""

    def draw(ax, positions):
        ax.fill_between(positions, 0, tn5_fwd, color="#A23B72", alpha=0.6, label="Forward")
        ax.fill_between(positions, 0, -tn5_rev, color="#2E86AB", alpha=0.8, label="Reverse")
        ax.axhline(0, color="black", linewidth=0.8)
        ax.legend(loc="upper right", frameon=False, fontsize=8)

    return _track_plot("Tn5 cut sites (n)", draw, len(tn5_fwd))


def create_read_start_sites_plot(tn5_fwd, tn5_rev, coverage_total):
    """Read start sites per position, from the recorded cut sites."""
    starts = tn5_fwd + tn5_rev
    label = "Read start sites (n)"
    if not starts.any():
        # --no-tn5 leaves no start record; fall back to depth and say so.
        starts = coverage_total
        label = "Total depth (n)"

    def draw(ax, positions):
        ax.fill_between(positions, 0, starts, color="#2E86AB", alpha=0.8)

    return _track_plot(label, draw, len(starts))


def create_tn5_insertion_context_plot(tn5_fwd, tn5_rev, reference):
    """Dinucleotide context of Tn5 insertion sites."""
    total_tn5 = tn5_fwd + tn5_rev
    dinuc_counts = {first + second: 0 for first in "ACGT" for second in "ACGT"}
    for pos in np.flatnonzero(total_tn5[:-1]):
        dinuc = reference[pos] + reference[pos + 1]
        if dinuc in dinuc_counts:
            dinuc_counts[dinuc] += int(total_tn5[pos])

    dinucs = sorted(dinuc_counts)
    counts = [dinuc_counts[d] for d in dinucs]
    total_cuts = sum(counts)
    if total_cuts == 0:
        return _empty_plot("No Tn5 cut data available", figsize=(8, 4))

    percentages = [count / total_cuts * 100 for count in counts]
    # Magenta for AT-only, blue for GC-only, purple for mixed contexts.
    palette = {0: "#A23B72", 1: "#9B59B6", 2: "#2E86AB"}
    colors = [palette[sum(base in "GC" for base in dinuc)] for dinuc in dinucs]

    fig, ax = plt.subplots(figsize=(8, 4))
    bars = ax.bar(dinucs, percentages, color=colors, alpha=0.8, edgecolor="black", linewidth=0.5)
    ax.set_ylim(0, max(percentages) * 1.1)
    _tidy(ax, "Dinucleotide context", "Tn5 insertion frequency (%)")
    ax.tick_params(labelsize=9)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    for bar, percentage in zip(bars, percentages, strict=True):
        if percentage > 0.5:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                percentage,
                f"{percentage:.1f}%",
                ha="center",
                va="bottom",
                fontsize=7,
            )
    return plot_to_base64(fig)


def create_depth_vs_coverage_plot(mean_depth, genome_coverage):
    """Coverage breadth against mean mtDNA depth."""
    mask = (mean_depth > 0) & (genome_coverage > 0)
    if not mask.any():
        return _empty_plot("No data available")

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(mean_depth[mask], genome_coverage[mask], s=20, color="black", edgecolors="none")
    ax.set_ylim(0, 1.05)
    _tidy(ax, "Mean mtDNA depth", "Coverage breadth")
    return plot_to_base64(fig)


def _log_scatter(x_values, y_values, x_label):
    """Log-log per-cell scatter against mtDNA depth."""
    mask = (y_values > 0) & (x_values > 0)
    if not mask.any():
        return _empty_plot("No data available\n(all values are zero)")

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(x_values[mask], y_values[mask], s=10, color="black", edgecolors="none")
    ax.set_xscale("log")
    ax.set_yscale("log")
    _tidy(ax, x_label, "mtDNA depth (log10)")
    ax.yaxis.set_major_locator(FixedLocator(_DEPTH_TICKS))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{int(value):,}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{int(value):,}"))
    return plot_to_base64(fig)


def create_depth_vs_fragments_plot(mean_depth, total_fragments):
    """mtDNA depth against total fragments per cell."""
    if total_fragments is None:
        return _empty_plot("No fragment counts available")
    return _log_scatter(total_fragments, mean_depth, "Total fragments (log10)")


def create_reads_vs_depth_plot(mean_depth, total_bases):
    """Estimated read count against mtDNA depth, for scRNA-seq."""
    # total_bases does not record read count; 150 bp is a useful scRNA-seq estimate.
    return _log_scatter(total_bases / 150.0, mean_depth, "Number of reads (log10)")


def _load_run_config(output_dir: Path) -> dict:
    """Read the run configuration the pipeline wrote beside the QC tables."""
    config_file = output_dir / "qc" / "run_config.json"
    if not config_file.exists():
        return {}
    try:
        with open(config_file) as f:
            return json.load(f)
    except (OSError, ValueError) as exc:
        logger.warning("Could not read %s: %s", config_file, exc)
        return {}


def _render_html(title, output_dir, run_config, data, sections):
    """Assemble the report page from already-rendered plot sections."""
    title = html.escape(title)
    parameters = run_config.get("parameters", {})
    parameter_rows = "\n".join(
        f"<tr><td>{name}</td><td>{value}</td></tr>" for name, value in sorted(parameters.items())
    )
    plot_html = "\n".join(sections)

    return f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title} - mgatk2 Report</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
            color: black;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: black;
            font-weight: bold;
            border-bottom: 3px solid black;
            padding-bottom: 10px;
            margin-bottom: 5px;
        }}
        .subtitle {{
            color: #666;
            font-size: 1.1em;
            font-style: italic;
            margin-bottom: 20px;
        }}
        h2 {{
            color: black;
            font-weight: bold;
            margin-top: 30px;
            border-left: 4px solid black;
            padding-left: 10px;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .stat-box {{
            background-color: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            border-left: 3px solid #2E86AB;
        }}
        .stat-label {{
            font-size: 0.85em;
            color: black;
            text-transform: uppercase;
        }}
        .stat-value {{
            font-size: 1.5em;
            font-weight: bold;
            color: black;
            margin-top: 5px;
        }}
        .plot {{
            margin: 20px 0;
            text-align: center;
        }}
        .plot img {{
            max-width: 100%;
            height: auto;
        }}
        .plot-grid {{
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 20px;
            margin: 20px 0;
        }}
        .plot-grid .plot {{
            margin: 0;
        }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            font-size: 0.9em;
            color: black;
            text-align: center;
        }}
        .version-info {{
            margin-top: 30px;
            padding: 15px;
            background-color: #f8f9fa;
            border-radius: 5px;
            font-size: 0.85em;
        }}
        .version-info h3 {{
            font-size: 1em;
            margin-top: 0;
            margin-bottom: 10px;
            border-left: none;
        }}
        .version-info table {{
            border-collapse: collapse;
            width: 100%;
        }}
        .version-info td {{
            padding: 3px 8px 3px 0;
            border-bottom: 1px solid #eee;
        }}
        .version-info td:first-child {{
            color: #666;
            width: 40%;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>
        <div class="subtitle">mgatk2 output analysis</div>
        <p><strong>Date:</strong> {datetime.now().strftime("%d/%m/%Y, %H:%M")}</p>
        <p><strong>Working directory:</strong> {os.getcwd()}</p>

        <h2>Summary statistics</h2>
        <div class="summary-grid">
            <div class="stat-box">
                <div class="stat-label">Total cells</div>
                <div class="stat-value">{data["n_cells"]:,}</div>
            </div>
            <div class="stat-box">
                <div class="stat-label">Cells passing QC</div>
                <div class="stat-value">{run_config.get("cells_passed_qc", "N/A")}</div>
            </div>
            <div class="stat-box">
                <div class="stat-label">Mean depth</div>
                <div class="stat-value">{data["mean_depth"].mean():.1f}×</div>
            </div>
        </div>

{plot_html}

        <div class="footer">
            Generated by mgatk2 v{run_config.get("mgatk_version", "unknown")} |
            Reference: {run_config.get("reference", "chrM")} |
            Output: {output_dir.name}
        </div>

        <div class="version-info">
            <h3>Run configuration</h3>
            <p>Python {sys.version.split()[0]} | input {run_config.get("input_bam", "unknown")}</p>
            <details>
                <summary style="cursor: pointer; font-weight: bold;">Parameters</summary>
                <table>{parameter_rows}</table>
            </details>
        </div>
    </div>
</body>
</html>
"""


def _plot_section(heading, plot, caption=""):
    caption_html = (
        f'\n            <p style="color: #666; font-size: 0.9em; margin-top: 10px;">{caption}</p>'
        if caption
        else ""
    )
    return f"""        <h2>{heading}</h2>
        <div class="plot">
            <img src="{plot}" alt="{heading}">{caption_html}
        </div>
"""


def _plot_pair(first_heading, first_plot, second_heading, second_plot):
    return f"""        <div class="plot-grid">
            <div>
                <h2>{first_heading}</h2>
                <div class="plot"><img src="{first_plot}" alt="{first_heading}"></div>
            </div>
            <div>
                <h2>{second_heading}</h2>
                <div class="plot"><img src="{second_plot}" alt="{second_heading}"></div>
            </div>
        </div>
"""


def generate_html_report(output_dir: Path, title: str, tn5: bool = True) -> Path:
    """Write the QC report: Tn5 transposition plots for ATAC, read starts otherwise."""
    output_dir = Path(output_dir)
    data = load_report_data(
        output_dir / "output" / "counts.h5", output_dir / "output" / "metadata.h5"
    )

    sections = [_plot_section("chrM coverage", create_coverage_plot(data["coverage_mean"]))]
    if tn5:
        sections += [
            _plot_section(
                "Tn5 transposition frequency",
                create_transposition_frequency_plot(data["tn5_fwd"], data["tn5_rev"]),
            ),
            _plot_section(
                "Tn5 insertion sequence context",
                create_tn5_insertion_context_plot(
                    data["tn5_fwd"], data["tn5_rev"], data["reference"]
                ),
                caption="magenta = AT-rich, blue = GC-rich, purple = mixed",
            ),
            _plot_pair(
                "Depth per cell",
                create_depth_vs_fragments_plot(data["mean_depth"], data["total_fragments"]),
                "chrM coverage",
                create_depth_vs_coverage_plot(data["mean_depth"], data["genome_coverage"]),
            ),
        ]
    else:
        sections += [
            _plot_section(
                "Read start sites",
                create_read_start_sites_plot(
                    data["tn5_fwd"], data["tn5_rev"], data["coverage_total"]
                ),
            ),
            _plot_pair(
                "Number of reads",
                create_reads_vs_depth_plot(data["mean_depth"], data["total_bases"]),
                "chrM coverage",
                create_depth_vs_coverage_plot(data["mean_depth"], data["genome_coverage"]),
            ),
        ]

    report_file = output_dir / "mgatk2_report.html"
    report_file.write_text(
        _render_html(title, output_dir, _load_run_config(output_dir), data, sections)
    )
    return report_file

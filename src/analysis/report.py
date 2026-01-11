"""HTML report generation"""

import base64
import logging
from io import BytesIO
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")

import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def plot_to_base64(fig):
    """Convert matplotlib figure to PNG"""
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode("utf-8")
    buf.close()
    plt.close(fig)
    return f"data:image/png;base64,{img_base64}"


def create_coverage_plot(coverage, positions=None):
    """Create coverage plot with log10 scale"""
    if coverage.ndim == 2:
        mean_coverage = coverage.mean(axis=1)
    else:
        mean_coverage = coverage

    if positions is None:
        positions = np.arange(1, len(mean_coverage) + 1)

    fig, ax = plt.subplots(figsize=(10, 3.5))

    ax.plot(positions, mean_coverage, linewidth=0.8, color="#2E86AB", alpha=1, zorder=2)

    ax.set_xlabel("chrM (bp)", fontsize=10, color="black")
    ax.set_ylabel("Mean depth", fontsize=10, color="black")
    ax.set_xlim(0, positions[-1])
    ax.tick_params(colors="black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")

    return plot_to_base64(fig)


def create_transposition_frequency_plot(counts_file):
    """Create mirrored transposition frequency plot"""
    with h5py.File(counts_file, "r") as f:
        # Sum Tn5 cut sites across all cells for each strand
        tn5_fwd = f["tn5_cuts_fwd"][:, :].sum(axis=1).astype(np.int64)
        tn5_rev = f["tn5_cuts_rev"][:, :].sum(axis=1).astype(np.int64)
        n_positions = len(tn5_fwd)

    positions = np.arange(1, n_positions + 1)

    fig, ax = plt.subplots(figsize=(10, 3.5))

    # Plot forward strand above x-axis (positive)
    ax.fill_between(
        positions,
        0,
        tn5_fwd,
        linewidth=0.5,
        color="#A23B72",
        alpha=0.6,
        label="Forward",
    )
    # Plot reverse strand below x-axis (negative)
    ax.fill_between(
        positions,
        0,
        -tn5_rev,
        linewidth=0.5,
        color="#2E86AB",
        alpha=0.8,
        label="Reverse",
    )

    ax.set_xlabel("chrM (bp)", fontsize=10, color="black")
    ax.set_ylabel("Tn5 cut sites (n)", fontsize=10, color="black")
    ax.set_xlim(0, positions[-1])
    ax.axhline(0, color="black", linewidth=0.8, linestyle="-")
    ax.tick_params(colors="black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")
    ax.legend(loc="upper right", frameon=False, fontsize=8)

    return plot_to_base64(fig)


def create_depth_vs_coverage_plot(metadata_file):
    """Coverage v s mtDNA depth plot"""
    with h5py.File(metadata_file, "r") as f:
        mean_depth = f["mean_depth"][:]
        genome_coverage = f["genome_coverage"][:]

    # Filter out zeros
    mask = (mean_depth > 0) & (genome_coverage > 0)
    mean_depth = mean_depth[mask]
    genome_coverage = genome_coverage[mask]

    if len(mean_depth) == 0:
        logger.warning("No valid data for depth vs coverage plot")
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.text(
            0.5,
            0.5,
            "No data available",
            ha="centre",
            va="centre",
            fontsize=14,
            color="gray",
        )
        ax.axis("off")
        return plot_to_base64(fig)

    fig, ax = plt.subplots(figsize=(6, 5))

    ax.scatter(mean_depth, genome_coverage, s=20, alpha=1, color="black", edgecolors="none")

    ax.set_ylabel("Coverage breadth (%)", fontsize=10, color="black")
    ax.set_xlabel("Mean mtDNA depth", fontsize=10, color="black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0, 105)

    return plot_to_base64(fig)
    return plot_to_base64(fig)


def create_depth_vs_fragments_plot(metadata_file):
    """depth vs total fragments plot"""
    with h5py.File(metadata_file, "r") as f:
        mtdna_depth = f["mean_depth"][:]
        total_fragments = f["barcode_metadata"]["total"][:]

    # Filter out zeros
    mask = (mtdna_depth > 0) & (total_fragments > 0)
    mtdna_depth = mtdna_depth[mask]
    total_fragments = total_fragments[mask]

    # Check if we have data
    if len(mtdna_depth) == 0:
        logger.warning("No valid data for depth vs fragments plot")
        # Return placeholder
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.text(
            0.5,
            0.5,
            "No data available\n(all values are zero)",
            ha="centre",
            va="centre",
            fontsize=14,
            color="gray",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        return plot_to_base64(fig)

    fig, ax = plt.subplots(figsize=(6, 5))

    ax.scatter(total_fragments, mtdna_depth, s=10, alpha=1, color="black", edgecolors="none")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Total fragments (log10)", fontsize=10, color="black")
    ax.set_ylabel("mtDNA depth (log10)", fontsize=10, color="black")
    ax.tick_params(colors="black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")

    # Set specific y-axis tick positions and labels
    from matplotlib.ticker import FixedLocator, FuncFormatter

    yticks = [1, 5, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500, 1000, 5000, 10000]
    ax.yaxis.set_major_locator(FixedLocator(yticks))

    # Format tick labels
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, p: f"{int(x):,}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, p: f"{int(y):,}"))

    return plot_to_base64(fig)


def generate_html_report(
    output_dir: Path,
    sample_name: str = "mgatk2",
    title: str | None = None,
    subtitle: str | None = None,
    working_directory: str | None = None,
    input_dir: str | None = None,
):
    """Generate HTML report with QC plots"""
    output_dir = Path(output_dir)

    # Auto-detect title from 10X run folder if not provided
    if title is None and input_dir is not None:
        input_path = Path(input_dir)
        # Check if this is a 10X structure (has 'outs' folder)
        if input_path.name == "outs":
            # Use parent folder name (the actual 10X run name)
            title = input_path.parent.name
        elif (input_path / "outs").exists():
            # Input dir is the 10X run folder itself
            title = input_path.name
        else:
            # Fallback to input directory name
            title = input_path.name

    # Auto-detect from working directory if input_dir not provided
    elif title is None and working_directory is not None:
        work_path = Path(working_directory)
        # Look for 10X structure indicators
        if work_path.name == "outs":
            title = work_path.parent.name
        elif (work_path / "outs").exists():
            title = work_path.name
        elif work_path.parent.name == "outs":
            title = work_path.parent.parent.name
        else:
            # Use working directory name
            title = work_path.name

    # Fallback to sample_name if still no title
    if title is None:
        title = sample_name

    counts_file = output_dir / "output" / "counts.h5"
    metadata_file = output_dir / "output" / "metadata.h5"

    if not counts_file.exists() or not metadata_file.exists():
        logger.error("Output files not found in %s", output_dir)
        return None

    # Load data and create plots
    with h5py.File(metadata_file, "r") as f:
        coverage = f["coverage"][:]
    coverage_plot = create_coverage_plot(coverage)

    transposition_plot = create_transposition_frequency_plot(counts_file)

    depth_plot = create_depth_vs_fragments_plot(metadata_file)

    depth_vs_coverage_plot = create_depth_vs_coverage_plot(metadata_file)

    # Load summary statistics
    with h5py.File(metadata_file, "r") as f:
        n_cells = len(f["mean_depth"][:])
        mean_depth_vals = f["mean_depth"][:]

        mean_depth = mean_depth_vals.mean()

    # Read summary file if exists
    summary_file = output_dir / "qc" / "summary.txt"
    summary_stats = {}
    if summary_file.exists():
        with open(summary_file) as f:
            for line in f:
                if ":" in line:
                    key, value = line.split(":", 1)
                    summary_stats[key.strip()] = value.strip()

    # Use current date/time in dd/mm/yyyy, HH:MM format
    from datetime import datetime

    run_date = datetime.now().strftime("%d/%m/%Y, %H:%M")

    # Set default title and subtitle if not provided
    if title is None:
        title = sample_name
    if subtitle is None:
        subtitle = "mgatk2 output analysis"

    # Generate HTML
    html_content = f"""
<!DOCTYPE html>
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
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>
        <div class="subtitle">{subtitle}</div>
        <p><strong>Date:</strong> {run_date}</p>
        {
        f"<p><strong>Working directory:</strong> {working_directory}</p>"
        if working_directory
        else ""
    }

        <h2>Summary statistics</h2>
        <div class="summary-grid">
            <div class="stat-box">
                <div class="stat-label">Total cells</div>
                <div class="stat-value">{n_cells:,}</div>
            </div>
            <div class="stat-box">
                <div class="stat-label">Cells passing QC</div>
                <div class="stat-value">{summary_stats.get("cells_passed_qc", "N/A")}</div>
            </div>
            <div class="stat-box">
                <div class="stat-label">Mean depth</div>
                <div class="stat-value">{mean_depth:.1f}×</div>
            </div>
        </div>

        <h2>chrM coverage</h2>
        <div class="plot">
            <img src="{coverage_plot}" alt="Coverage Plot">
        </div>

        <h2>Tn5 transposition frequency</h2>
        <div class="plot">
            <img src="{transposition_plot}" alt="Transposition Frequency Plot">
        </div>

        <div class="plot-grid">
            <div>
                <h2>Depth per cell</h2>
                <div class="plot">
                    <img src="{depth_plot}" alt="Depth vs Fragments Plot">
                </div>
            </div>
            <div>
                <h2>chrM coverage</h2>
                <div class="plot">
                    <img src="{depth_vs_coverage_plot}" alt="mtDNA Depth vs Genome Coverage">
                </div>
            </div>
        </div>

        <div class="footer">
            Generated by mgatk2 v{summary_stats.get("mgatk_version", "1.0.0")} |
            Reference: {summary_stats.get("reference", "chrM")} |
            Output: {output_dir.name}
        </div>
    </div>
</body>
</html>
"""

    # Write HTML file
    report_file = output_dir / "mgatk2_report.html"
    with open(report_file, "w") as f:
        f.write(html_content)

    return report_file

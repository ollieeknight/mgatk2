"""HTML report generation"""

import base64
import logging
import subprocess
from io import BytesIO
from pathlib import Path

import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")

logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)

logger = logging.getLogger(__name__)


def get_environment_info():
    """Collect conda and pip package versions for the report"""
    env_info = {"conda": "", "pip": ""}

    # Try to get conda list
    try:
        result = subprocess.run(["conda", "list"], capture_output=True, text=True, timeout=30)
        if result.returncode == 0:
            env_info["conda"] = result.stdout
    except (subprocess.TimeoutExpired, FileNotFoundError, Exception) as e:
        logger.debug("Could not get conda list: %s", e)
        env_info["conda"] = "conda not available"

    # Try to get pip list
    try:
        result = subprocess.run(["pip", "list"], capture_output=True, text=True, timeout=30)
        if result.returncode == 0:
            env_info["pip"] = result.stdout
    except (subprocess.TimeoutExpired, FileNotFoundError, Exception) as e:
        logger.debug("Could not get pip list: %s", e)
        env_info["pip"] = "pip not available"

    return env_info


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


def create_read_start_sites_plot(metadata_file):
    """Create read start sites plot for scRNA-seq data (positive y-axis only)"""
    with h5py.File(metadata_file, "r") as f:
        # Sum coverage across all cells to get total read starts at each position
        coverage = f["coverage"][:, :].sum(axis=1).astype(np.int64)
        n_positions = len(coverage)

    positions = np.arange(1, n_positions + 1)

    fig, ax = plt.subplots(figsize=(10, 3.5))

    # Plot read start sites as positive values only
    ax.fill_between(
        positions,
        0,
        coverage,
        linewidth=0.5,
        color="#2E86AB",
        alpha=0.8,
    )

    ax.set_xlabel("chrM (bp)", fontsize=10, color="black")
    ax.set_ylabel("Read start sites (n)", fontsize=10, color="black")
    ax.set_xlim(0, positions[-1])
    ax.tick_params(colors="black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")

    return plot_to_base64(fig)


def create_tn5_insertion_context_plot(counts_file, metadata_file):
    """Create bar plot showing dinucleotide context of Tn5 insertion sites"""
    with h5py.File(counts_file, "r") as f:
        tn5_fwd = f["tn5_cuts_fwd"][:, :].sum(axis=1).astype(np.int64)
        tn5_rev = f["tn5_cuts_rev"][:, :].sum(axis=1).astype(np.int64)

    with h5py.File(metadata_file, "r") as f:
        refallele = f["reference"][:]
        if isinstance(refallele[0], bytes):
            refallele = [x.decode() for x in refallele]
        else:
            refallele = list(refallele)

    # Combine forward and reverse Tn5 cuts
    total_tn5 = tn5_fwd + tn5_rev

    # Count dinucleotide contexts at Tn5 insertion sites
    dinuc_counts = {}
    bases = ["A", "C", "G", "T"]
    for b1 in bases:
        for b2 in bases:
            dinuc_counts[f"{b1}{b2}"] = 0

    # For each position with Tn5 cuts, get the dinucleotide context
    for pos in range(len(total_tn5) - 1):
        if total_tn5[pos] > 0:
            base1 = refallele[pos]
            base2 = refallele[pos + 1]
            dinuc = f"{base1}{base2}"
            if dinuc in dinuc_counts:
                dinuc_counts[dinuc] += total_tn5[pos]

    # Sort by dinucleotide for consistent ordering
    dinucs = sorted(dinuc_counts.keys())
    counts = [dinuc_counts[d] for d in dinucs]

    # Calculate percentages
    total_cuts = sum(counts)
    if total_cuts == 0:
        logger.warning("No Tn5 cuts found for insertion context plot")
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(
            0.5,
            0.5,
            "No Tn5 cut data available",
            ha="center",
            va="center",
            fontsize=14,
            color="gray",
        )
        ax.axis("off")
        return plot_to_base64(fig)

    percentages = [(c / total_cuts) * 100 for c in counts]

    # Create bar plot
    fig, ax = plt.subplots(figsize=(8, 4))

    # Color by base composition
    colors = []
    for dinuc in dinucs:
        gc_content = (dinuc.count("G") + dinuc.count("C")) / 2
        if gc_content == 0:
            colors.append("#A23B72")  # AT-rich (magenta)
        elif gc_content == 1:
            colors.append("#2E86AB")  # GC-rich (blue)
        else:
            colors.append("#9B59B6")  # Mixed (purple)

    bars = ax.bar(dinucs, percentages, color=colors, alpha=0.8, edgecolor="black", linewidth=0.5)

    ax.set_xlabel("Dinucleotide context", fontsize=10, color="black")
    ax.set_ylabel("Tn5 insertion frequency (%)", fontsize=10, color="black")
    ax.set_ylim(0, max(percentages) * 1.1)
    ax.tick_params(colors="black", labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")
    plt.xticks(rotation=45, ha="right")

    # Add value labels on bars
    for bar, pct in zip(bars, percentages, strict=True):
        height = bar.get_height()
        if height > 0.5:  # Only show label if bar is visible
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                height,
                f"{pct:.1f}%",
                ha="center",
                va="bottom",
                fontsize=7,
                color="black",
            )

    plt.tight_layout()

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


def create_reads_vs_depth_plot(metadata_file):
    """Number of reads vs mtDNA depth plot for scRNA-seq"""
    with h5py.File(metadata_file, "r") as f:
        mtdna_depth = f["mean_depth"][:]
        total_bases = f["total_bases"][:]

    # Estimate number of reads from total_bases (assuming ~150bp read length for scRNA-seq)
    n_reads = total_bases / 150.0

    # Filter out zeros
    mask = (mtdna_depth > 0) & (n_reads > 0)
    mtdna_depth = mtdna_depth[mask]
    n_reads = n_reads[mask]

    # Check if we have data
    if len(mtdna_depth) == 0:
        logger.warning("No valid data for reads vs depth plot")
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

    ax.scatter(n_reads, mtdna_depth, s=10, alpha=1, color="black", edgecolors="none")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of reads (log10)", fontsize=10, color="black")
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

    tn5_context_plot = create_tn5_insertion_context_plot(counts_file, metadata_file)

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

    # Get environment info for version section
    env_info = get_environment_info()

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
        .version-info {{
            margin-top: 30px;
            padding: 15px;
            background-color: #f8f9fa;
            border-radius: 5px;
            font-size: 0.8em;
        }}
        .version-info h3 {{
            font-size: 1em;
            margin-top: 0;
            margin-bottom: 10px;
            border-left: none;
        }}
        .version-info pre {{
            background-color: white;
            padding: 10px;
            border-radius: 3px;
            overflow-x: auto;
            white-space: pre-wrap;
            word-wrap: break-word;
            max-height: 300px;
            overflow-y: auto;
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

        <h2>Tn5 insertion sequence context</h2>
        <div class="plot">
            <img src="{tn5_context_plot}" alt="Tn5 Insertion Context Plot">
            <p style="color: #666; font-size: 0.9em; margin-top: 10px;">
                magenta = AT-rich, blue = GC-rich, purple = mixed
            </p>
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

        <div class="version-info">
            <h3>Environment Information</h3>
            <details>
                <summary style="cursor: pointer; font-weight: bold;">Conda Packages</summary>
                <pre>{env_info["conda"]}</pre>
            </details>
            <details>
                <summary style="cursor: pointer; font-weight: bold;">Pip Packages</summary>
                <pre>{env_info["pip"]}</pre>
            </details>
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


def generate_scrna_html_report(
    output_dir: Path,
    sample_name: str = "mgatk2",
    title: str | None = None,
    subtitle: str | None = None,
    working_directory: str | None = None,
    input_dir: str | None = None,
):
    """Generate HTML report for scRNA-seq data (without singlecell.csv)"""
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

    # Use read start sites instead of Tn5 transposition
    read_starts_plot = create_read_start_sites_plot(metadata_file)

    # Use reads vs depth instead of fragments vs depth
    reads_depth_plot = create_reads_vs_depth_plot(metadata_file)

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
        subtitle = "mgatk2 scRNA-seq output analysis"

    # Get environment info for version section
    env_info = get_environment_info()

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
        .version-info {{
            margin-top: 30px;
            padding: 15px;
            background-color: #f8f9fa;
            border-radius: 5px;
            font-size: 0.8em;
        }}
        .version-info h3 {{
            font-size: 1em;
            margin-top: 0;
            margin-bottom: 10px;
            border-left: none;
        }}
        .version-info pre {{
            background-color: white;
            padding: 10px;
            border-radius: 3px;
            overflow-x: auto;
            white-space: pre-wrap;
            word-wrap: break-word;
            max-height: 300px;
            overflow-y: auto;
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

        <h2>Read start sites</h2>
        <div class="plot">
            <img src="{read_starts_plot}" alt="Read Start Sites Plot">
        </div>

        <div class="plot-grid">
            <div>
                <h2>Number of reads</h2>
                <div class="plot">
                    <img src="{reads_depth_plot}" alt="Reads vs mtDNA Depth">
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

        <div class="version-info">
            <h3>Environment Information</h3>
            <details>
                <summary style="cursor: pointer; font-weight: bold;">Conda Packages</summary>
                <pre>{env_info["conda"]}</pre>
            </details>
            <details>
                <summary style="cursor: pointer; font-weight: bold;">Pip Packages</summary>
                <pre>{env_info["pip"]}</pre>
            </details>
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

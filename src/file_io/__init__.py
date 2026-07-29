"""IO module for mgatk2"""

from .barcode_extraction import extract_barcodes_from_bam
from .formats import write_cell_stats, write_run_summary
from .writers import IncrementalHDF5Writer, IncrementalTextWriter

__all__ = [
    "IncrementalHDF5Writer",
    "IncrementalTextWriter",
    "write_cell_stats",
    "write_run_summary",
    "extract_barcodes_from_bam",
]

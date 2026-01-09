"""IO module for mgatk2"""

from .barcode_extraction import extract_barcodes_from_bam
from .formats import (
    detect_delimiter,
    get_compression_type,
    validate_file_format,
    write_cell_stats,
    write_position_stats,
    write_reference_fasta,
    write_run_summary,
    write_tsv_file,
)
from .main import save_mgatk_outputs
from .writers import IncrementalHDF5Writer, IncrementalTextWriter, write_parameters_json

__all__ = [
    "save_mgatk_outputs",
    "IncrementalHDF5Writer",
    "IncrementalTextWriter",
    "write_parameters_json",
    "write_cell_stats",
    "write_position_stats",
    "write_run_summary",
    "write_reference_fasta",
    "write_tsv_file",
    "validate_file_format",
    "detect_delimiter",
    "get_compression_type",
    "extract_barcodes_from_bam",
]

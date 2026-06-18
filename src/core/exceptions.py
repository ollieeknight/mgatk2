"""Exceptions for mgatk2."""


class MgatkError(Exception):
    """Base exception for mgatk2 errors."""


class InvalidInputError(MgatkError):
    """Raised when input files are invalid or missing"""


class ProcessingError(MgatkError):
    """Raised when pipeline processing fails"""


class BAMReadError(ProcessingError):
    """Raised when BAM file reading fails"""

    def __init__(self, bam_path: str, message: str):
        self.bam_path = bam_path
        super().__init__(f"BAM read error for {bam_path}: {message}")


class NoChrMReadsError(BAMReadError):
    """Raised when BAM file contains no mitochondrial reads"""

    def __init__(self, bam_path: str, available_chromosomes: list[str]):
        self.available_chromosomes = available_chromosomes
        chrs = ", ".join(available_chromosomes[:10])
        super().__init__(
            bam_path,
            f"No mitochondrial chromosome (chrM, MT, or M) found.\n"
            f"Available: {chrs}{'...' if len(available_chromosomes) > 10 else ''}",
        )


class NoBarcodeTagsError(BAMReadError):
    """Raised when BAM file contains no reads with the specified barcode tag"""

    def __init__(self, bam_path: str, barcode_tag: str, total_reads_checked: int):
        self.barcode_tag = barcode_tag
        self.total_reads_checked = total_reads_checked
        super().__init__(
            bam_path,
            (
                f"No reads with barcode tag '{barcode_tag}' found "
                f"(checked {total_reads_checked:,} reads).\n"
                f"This may not be a single-cell BAM file, or wrong tag specified."
            ),
        )


class BAMFormatError(InvalidInputError):
    """Raised when BAM file is corrupted or invalid format"""

    def __init__(self, bam_path: str, details: str = ""):
        message = f"BAM file appears corrupted or is not a valid BAM format: {bam_path}"
        if details:
            message += f"\n{details}"
        super().__init__(message)

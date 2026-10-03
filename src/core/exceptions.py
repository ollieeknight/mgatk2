"""Exceptions for mgatk2."""


class MgatkError(Exception):
    """Base exception for mgatk2 errors."""


class InvalidInputError(MgatkError):
    """Raised when input files are invalid or missing"""


class ProcessingError(MgatkError):
    """Raised when pipeline processing fails"""

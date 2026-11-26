"""turnpike: reconstruct integer points on a line from their pairwise distances.

    >>> from turnpike import solve
    >>> rep = solve([1, 4, 10, 12, 17, 3, 9, 11, 16, 6, 8, 13, 2, 7, 5])   # any order
    >>> rep.status, rep.points

See README.md for the sections and their guarantees.
"""
from .distances import Instance, InvalidInput, canonical, distances_of, verify
from .options import Options
from .report import Report
from .solver import solve

__all__ = ["solve", "Options", "Report", "Instance", "InvalidInput", "verify", "canonical", "distances_of"]
__version__ = "0.1.0"

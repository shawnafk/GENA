"""Source access APIs for external observation products."""

from .indices import read_indices, parse_indices
from .spedas import parse_spedas, read_spedas

__all__ = ["parse_spedas", "read_spedas", "parse_indices", "read_indices"]

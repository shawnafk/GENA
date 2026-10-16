"""Public source loaders, parsers, and time-series assembly."""

from load.load import (
    get_goes_loc,
    get_mms_loc,
    get_queqiao_loc,
    get_themis_loc,
    load_gena,
    load_indices,
    load_spedas,
    load,
    load_pos,
)
from load.api.indices import parse_indices, read_indices
from load.api.spedas import parse_spedas, read_spedas

__all__ = [
    "get_themis_loc",
    "get_mms_loc",
    "get_goes_loc",
    "get_queqiao_loc",
    "load",
    "load_pos",
    "load_gena",
    "load_indices",
    "load_spedas",
    "parse_indices",
    "parse_spedas",
    "read_indices",
    "read_spedas",
]

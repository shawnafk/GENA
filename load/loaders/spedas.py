"""Shared SPEDAS loading pipeline and mission-specific scientific helpers."""

import logging
import os
import subprocess

import numpy as np
import pandas as pd

from common.config import get_configured_paths
from common.time import Time
from load.api.spedas import (
    OMNI_FIELDS, parse_spedas, read_spedas,
)
from load.api.public import _clip_variable_dict, _merge_variable_dicts, assemble
from common.frame import concat
from load.api.public import normalize_interval

logger = logging.getLogger(__name__)
OMNI_VARIABLES = list(OMNI_FIELDS["hro"]["default"])
_mms_configured = False


def load_spedas(
    mission, trange, product=None, *, variables=None, probe=None, level="l2",
    prefix=None, download=None, return_type="dataframe", variable_layouts=None,
    **kwargs,
):
    """Load SPEDAS observations through one interval/assembly pipeline."""
    if return_type not in {"dataframe", "dict"}:
        raise ValueError("return_type must be 'dataframe' or 'dict'")
    mission = str(mission).lower()
    if mission not in {"themis", "mms", "omni"}:
        raise ValueError("mission must be 'themis', 'mms', or 'omni'")
    if product is None and mission == "omni":
        product = "hro"
    if mission == "omni" and product != "hro":
        raise ValueError("OMNI loader currently supports product='hro' only.")
    if mission == "themis" and product == "esd":
        if level != "l2":
            raise ValueError("THEMIS ESD is available only at level='l2'")
        datatype = kwargs.setdefault("datatype", "peif")
        if not isinstance(datatype, str) or not datatype:
            raise ValueError("THEMIS ESD datatype must be a non-empty string")
    probe = ("a" if mission == "themis" else 1) if probe is None and mission != "omni" else probe
    chunks, raw_chunks = [], []
    interval = normalize_interval(trange)
    for interval in [interval]:
        dates = [value.isoformat() for value in interval]
        raw, interval_layouts = read_spedas(
            mission, dates, product=product, variables=variables, probe=probe,
            level=level, prefix=prefix, download=download,
            variable_layouts=variable_layouts, **kwargs,
        )
        raw = _clip_variable_dict(raw, interval, layouts=interval_layouts)
        if return_type == "dict":
            raw_chunks.append(raw)
            continue
        parsed = parse_spedas(raw, interval_layouts)
        provenance = {"source": mission, "product": product, "level": "hro" if mission == "omni" else level}
        if mission != "omni":
            provenance["probe"] = probe
        frame = assemble(parsed, **provenance)
        if mission == "themis" and product == "esd":
            frame.attrs["datatype"] = kwargs["datatype"]
        chunks.append(frame)
    if return_type == "dataframe":
        return concat(chunks)
    return _merge_variable_dicts(raw_chunks, layouts=variable_layouts)




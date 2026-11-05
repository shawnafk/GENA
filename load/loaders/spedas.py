"""Shared SPEDAS loading pipeline."""

import logging

from load.api.spedas import SPEDAS_OPTIONS, default_level, normalize_probe, parse_spedas, read_spedas
from load.api.variable import _clip_variable_dict, assemble
from common.time import Time

logger = logging.getLogger(__name__)


def load_spedas(
    mission, trange, product=None, *, variables=None, probe=None, level=None,
    prefix=None, download=None, return_type="dataframe", **kwargs,
):
    """Load one SPEDAS product over one inclusive UTC interval.

    ``variable_layouts`` is an internal backend override; it is accepted only
    through ``kwargs`` and is not part of the documented loader arguments.
    """
    variable_layouts = kwargs.pop("variable_layouts", None)
    if return_type not in {"dataframe", "dict"}:
        raise ValueError("return_type must be 'dataframe' or 'dict'")
    mission = str(mission).lower()
    if mission not in SPEDAS_OPTIONS:
        raise ValueError(f"mission must be one of {', '.join(SPEDAS_OPTIONS)}")
    products = SPEDAS_OPTIONS[mission]["products"]
    if product is None and len(products) == 1:
        product = products[0]
    if product is not None:
        product = str(product).lower()
    if mission == "omni" and product != "hro":
        raise ValueError("OMNI loader currently supports product='hro' only.")
    level = default_level(mission, product) if level is None else level
    if mission == "themis" and product == "esd":
        if level != "l2":
            raise ValueError("THEMIS ESD is available only at level='l2'")
        datatype = kwargs.setdefault("datatype", "peif")
        if not isinstance(datatype, str) or not datatype:
            raise ValueError("THEMIS ESD datatype must be a non-empty string")
    probe = normalize_probe(mission, probe)

    interval = Time(trange).datetime
    dates = [value.isoformat() for value in interval]
    raw, interval_layouts = read_spedas(
        mission, dates, product=product, variables=variables, probe=probe,
        level=level, prefix=prefix, download=download,
        variable_layouts=variable_layouts, **kwargs,
    )
    raw = _clip_variable_dict(raw, interval, layouts=interval_layouts)
    if return_type == "dict":
        return raw

    parsed = parse_spedas(raw, interval_layouts)
    provenance = {
        "source": mission,
        "product": product,
        "level": "hro" if mission == "omni" else level,
    }
    if mission != "omni":
        provenance["probe"] = probe
    frame = assemble(parsed, **provenance)
    if not parsed:
        reason = (
            f"{mission} {product} returned no variables for "
            f"{dates[0]}..{dates[1]}"
        )
        frame.attrs["empty_reason"] = reason
        logger.warning("load_spedas: %s", reason)
    if mission == "themis" and product == "esd":
        frame.attrs["datatype"] = kwargs["datatype"]
    return frame

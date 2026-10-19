"""One index-loading pipeline; source access and field rules live in the API."""

from load.api.indices import read_indices, parse_indices
from load.api.public import normalize_interval
from load.api.public import _merge_variable_dicts, assemble
from common.frame import concat

__all__ = ["load_indices"]


def load_indices(trange, product=None, *, prefix=None, download="local_first", return_type="dataframe", timeout=30):
    """Load Kyoto/GFZ indices over one or more inclusive UTC intervals.

    Defaults to AE/AL/AU. Pass a name or list of names for other products;
    see ``load.api.indices.INDICES_FIELDS``. Native sampling and units are
    retained; different cadences are outer-joined without interpolation.
    ``read_indices`` owns cache, download policy, field selection and layouts.
    Dictionary output keeps the existing single-product x/y form, or a
    name-to-record mapping for multiple products.
    """
    if return_type not in {"dataframe", "dict"}:
        raise ValueError("return_type must be 'dataframe' or 'dict'")
    chunks, raw_chunks = [], []
    interval = normalize_interval(trange)
    for interval in [interval]:
        data, layouts = read_indices(interval, product, prefix=prefix, download=download, timeout=timeout)
        if return_type == "dict":
            raw_chunks.append(data)
        else:
            names = list(data)
            chunks.append(assemble(
                parse_indices(data, layouts), source="indices",
                product=names[0] if len(names) == 1 else names,
            ))
    if return_type == "dataframe":
        return concat(chunks)
    records = _merge_variable_dicts(raw_chunks)
    if len(records) == 1:
        return next(iter(records.values()))
    return records

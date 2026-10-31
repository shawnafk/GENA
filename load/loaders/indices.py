"""One index-loading pipeline; source access and field rules live in the API."""

from load.api.indices import read_indices, parse_indices
from load.api.variable import assemble
from common.time import Time

__all__ = ["load_indices"]


def load_indices(trange, product=None, *, prefix=None, download="local_first", return_type="dataframe", timeout=30):
    """Load Kyoto/GFZ indices over one inclusive UTC interval.

    Defaults to AE/AL/AU. Pass a name or list of names for other products;
    see ``load.api.indices.INDICES_FIELDS``. Native sampling and units are
    retained; different cadences are outer-joined without interpolation.
    ``read_indices`` owns cache, download policy, field selection and layouts.
    Dictionary output keeps the existing single-product x/y form, or a
    name-to-record mapping for multiple products.
    """
    if return_type not in {"dataframe", "dict"}:
        raise ValueError("return_type must be 'dataframe' or 'dict'")
    interval = Time(trange).datetime
    data, layouts = read_indices(
        interval, product, prefix=prefix, download=download, timeout=timeout
    )
    if return_type == "dict":
        if len(data) == 1:
            return next(iter(data.values()))
        return data
    names = list(data)
    return assemble(
        parse_indices(data, layouts),
        source="indices",
        product=names[0] if len(names) == 1 else names,
    )

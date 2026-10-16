"""PySPEDAS adapter for Earth-centred Cartesian component transforms."""

import numpy as np
from common.time import Time

FRAMES = {"GSE", "GSEQ", "GSM", "SM", "GEI", "GEO", "MAG", "J2000", "HEE", "HAE", "HEEQ"}


def transform(values, times, *, source, target, unit, origin):
    """Call cotrans with arrays, never creating tplot variables."""
    if origin != 'earth':
        raise NotImplementedError("PySPEDAS adapter requires origin='earth'")
    source = {'ECI2000': 'J2000'}.get(source, source)
    target = {'ECI2000': 'J2000'}.get(target, target)
    if source not in FRAMES or target not in FRAMES:
        raise NotImplementedError(f"PySPEDAS does not support {source} → {target}")
    from pyspedas import cotrans

    if not len(values):
        return values.copy()
    result = cotrans(time_in=Time(times).to_timestamp(), data_in=values,
                     coord_in=source.lower(), coord_out=target.lower())
    if result is None or np.shape(result) != values.shape:
        raise ValueError(f"PySPEDAS failed to transform {source} → {target}")
    return np.asarray(result, dtype=float)

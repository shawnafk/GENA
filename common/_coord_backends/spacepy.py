"""SpacePy adapter for Earth-centred Cartesian component transforms."""

import numpy as np
from common.time import Time


def transform(values, times, *, source, target, unit, origin):
    """Preserve numeric units for Cartesian rotations; do not translate origins."""
    if origin != "earth":
        raise NotImplementedError("SpacePy standard frames require origin='earth'")
    from spacepy.coordinates import Coords, SYSAXES_TYPES

    for frame in (source, target):
        if frame not in SYSAXES_TYPES or SYSAXES_TYPES[frame]['car'] is None:
            raise NotImplementedError(f"SpacePy Cartesian frame {frame!r} is unsupported")
    if not len(values):
        return values.copy()
    coordinates = Coords(values, source, 'car', ticks=Time(times).to_tick(),
                         use_irbem=False)
    return np.asarray(coordinates.convert(target, 'car').data, dtype=float)

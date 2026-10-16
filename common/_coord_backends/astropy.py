"""Astropy native position transforms, with explicit units and origins.

GEO maps to ITRS. SpacePy ECI2000/J2000 is deliberately not aliased to
GCRS or ICRS; GSE/GSM require a separately validated frame extension.
"""

import numpy as np

FRAMES = {'ITRS': 'earth', 'GCRS': 'earth', 'CIRS': 'earth', 'TEME': 'earth',
          'ICRS': 'solar_system_barycenter', 'GALACTIC': 'solar_system_barycenter'}


def transform(values, times, *, source, target, unit, origin):
    """Transform native positions, rejecting origin changes and unitless vectors."""
    source = {'GEO': 'ITRS'}.get(source, source)
    target = {'GEO': 'ITRS'}.get(target, target)
    if source not in FRAMES or target not in FRAMES:
        raise NotImplementedError(
            f"Astropy native backend does not support {source} → {target}; "
            "GSE/GSM/ECI2000 are not native aliases")
    if FRAMES[source] != origin or FRAMES[target] != origin:
        raise NotImplementedError("Astropy adapter does not translate coordinate origins")
    if unit is None:
        raise ValueError("Astropy position transforms require an explicit length unit")
    import astropy.units as u
    from astropy import coordinates as ac
    from astropy.time import Time
    from astropy.utils import iers
    from common.constants import RE

    length_unit = u.Unit(float(RE) * u.m) if unit.lower() == 're' else u.Unit(unit)
    if not length_unit.is_equivalent(u.m):
        raise ValueError("Astropy position transforms require a length unit")
    if not len(values):
        return values.copy()
    date = Time(times, format='datetime64', scale='utc')
    classes = {'ITRS': ac.ITRS, 'GCRS': ac.GCRS, 'CIRS': ac.CIRS,
               'TEME': ac.TEME, 'ICRS': ac.ICRS, 'GALACTIC': ac.Galactic}
    def frame(name):
        return classes[name](**({'obstime': date} if FRAMES[name] == 'earth' else {}))
    representation = ac.CartesianRepresentation(values.T * length_unit)
    # Use packaged IERS tables, keeping a coordinate conversion free of network I/O.
    with iers.conf.set_temp('auto_download', False):
        result = frame(source).realize_frame(representation).transform_to(frame(target))
    return np.asarray(result.cartesian.xyz.to_value(length_unit).T, dtype=float)

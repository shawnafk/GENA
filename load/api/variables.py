"""Clip decoded backend variable dictionaries without changing their layout."""

from __future__ import annotations

from collections.abc import Mapping
import re

import numpy as np

from common.time import Time


def _clip_variable_dict(data, trange, *, layouts=None):
    """Clip decoded backend variables to an inclusive time interval."""
    start, end = Time(trange).datetime
    clipped = {}
    for name, variable in data.items():
        if not isinstance(variable, Mapping) or "x" not in variable or "y" not in variable:
            continue
        times = Time(variable["x"]).datetime
        mask = np.asarray((times >= start) & (times <= end))
        item = dict(variable)
        item["x"] = times[mask]
        layout = (layouts or {}).get(name, {})
        for key, source in variable.items():
            if key != "y" and not re.fullmatch(r"v\d*", key):
                continue
            values = np.asanyarray(source)
            time_axis = (
                layout.get("time_axis", 0)
                if key == "y"
                else layout.get("axes", {}).get(key, {}).get(
                    "time_axis", None if values.ndim <= 1 else 0
                )
            )
            if time_axis is None:
                item[key] = values
                continue
            if (
                not isinstance(time_axis, int)
                or not 0 <= time_axis < values.ndim
                or values.shape[time_axis] != len(times)
            ):
                raise ValueError(
                    f"Variable {name!r}, field {key!r}: shape {values.shape} does not match {len(times)} times on axis {time_axis}"
                )
            item[key] = np.take(values, np.flatnonzero(mask), axis=time_axis)
        clipped[name] = item
    return clipped

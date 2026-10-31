"""Assemble normalized variable descriptions into time-series frames."""

from __future__ import annotations

from copy import deepcopy

import numpy as np
import pandas as pd

from common.frame import TimeSeriesDataFrame
from common.time import Time


def assemble(variables, **provenance):
    """Assemble normalized variable descriptions into a time-series frame."""
    frames, metadata, used = [], {}, {"Datetime", "Timestamp"}
    for name, variable in variables.items():
        dates = Time(variable["time"]).datetime
        if dates.hasnans:
            raise ValueError(f"Variable {name!r} has invalid observation times")
        try:
            values = np.asarray(variable["values"])
        except ValueError as exc:
            raise ValueError(f"Variable {name!r} has irregular sample shapes") from exc
        time_axis = variable["time_axis"]
        if not isinstance(time_axis, int) or not 0 <= time_axis < values.ndim:
            raise ValueError(
                f"Variable {name!r}: invalid time axis {time_axis} for shape {values.shape}"
            )
        values = np.moveaxis(values, time_axis, 0)
        if len(values) != len(dates):
            raise ValueError(
                f"Variable {name!r}: shape {values.shape} does not match {len(dates)} times"
            )
        kind, names = variable["kind"], list(variable["columns"])
        if kind == "scalar" and values.ndim == 1 and len(names) == 1:
            columns = {names[0]: values}
        elif kind == "vector" and values.ndim == 2 and values.shape[1] == len(names):
            columns = {column: values[:, i] for i, column in enumerate(names)}
        elif kind == "array" and values.ndim >= 2 and len(names) == 1:
            cells = np.empty(len(values), dtype=object)
            cells[:] = [np.asarray(value).copy() for value in values]
            columns = {names[0]: cells}
        else:
            raise ValueError(
                f"Variable {name!r}: {kind!r} layout {names} is incompatible with shape {values.shape}"
            )
        info = {
            "kind": kind,
            "columns": names,
            "sample_shape": values.shape[1:],
            "units": variable.get("units"),
            "attrs": deepcopy(variable.get("attrs", {})),
            "axes": {},
        }
        for axis_name, spec in variable.get("axes", {}).items():
            axis = np.asarray(spec["values"])
            axis_time = spec["time_axis"]
            if axis_time is None:
                samples = np.broadcast_to(axis, (len(dates), *axis.shape))
            else:
                if not isinstance(axis_time, int) or not 0 <= axis_time < axis.ndim:
                    raise ValueError(
                        f"Variable {name!r}, axis {axis_name!r}: invalid time axis for shape {axis.shape}"
                    )
                samples = np.moveaxis(axis, axis_time, 0)
                if len(samples) != len(dates):
                    raise ValueError(
                        f"Variable {name!r}, axis {axis_name!r}: shape {axis.shape} does not match {len(dates)} times"
                    )
            column = spec["column"]
            if column in columns:
                raise ValueError(f"Variable {name!r}: duplicate column {column!r}")
            axis_cells = np.empty(len(samples), dtype=object)
            axis_cells[:] = [np.asarray(value).copy() for value in samples]
            columns[column] = axis_cells
            info["axes"][axis_name] = {
                "column": column,
                "static": axis_time is None,
                "sample_shape": samples.shape[1:],
            }
            if axis_time is None:
                info["axes"][axis_name]["values"] = axis.tolist()
            if spec.get("role") == "energy":
                info["energy_column"] = column
        if len(set(names)) != len(names) or used.intersection(columns):
            raise ValueError(f"Variable {name!r}: duplicate or reserved output columns")
        used.update(columns)
        frames.append(
            pd.DataFrame({"Datetime": dates, **columns})
            .drop_duplicates("Datetime", keep="first")
            .set_index("Datetime")
        )
        metadata[name] = info
    if frames:
        joined = pd.concat(frames, axis=1, join="outer", sort=False).sort_index()
        dates = Time(joined.index).datetime
        data = {"Datetime": dates, "Timestamp": dates.asi8.astype(float) / 1e9}
        data.update({name: joined[name].to_numpy() for name in joined})
    else:
        dates = Time([]).datetime
        data = {"Datetime": dates, "Timestamp": dates.asi8.astype(float) / 1e9}
    frame = TimeSeriesDataFrame(data)
    frame.attrs = {**provenance, "variables": metadata}
    return frame

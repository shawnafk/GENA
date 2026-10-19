"""Public time-series assembly and backend dictionary operations."""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
import re

import numpy as np
import pandas as pd

from common.frame import TimeSeriesDataFrame
from common.time import Time


def normalize_interval(trange):
    """Normalize one inclusive interval to two ordered UTC timestamps."""
    dates = Time(trange).datetime
    if isinstance(dates, pd.Timestamp) or len(dates) != 2:
        raise ValueError("trange must contain exactly two time values")
    if dates[0] > dates[1]:
        raise ValueError("trange start must not be later than end")
    return dates


def _time_sample_cells(values):
    cells = np.empty(len(values), dtype=object)
    cells[:] = [np.asarray(value).copy() for value in values]
    return cells


def _time_columns(times):
    dates = Time(times).datetime
    return {"Datetime": dates, "Timestamp": dates.asi8.astype(float) / 1e9}


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
            columns = {names[0]: _time_sample_cells(values)}
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
            columns[column] = _time_sample_cells(samples)
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
        data = _time_columns(joined.index)
        data.update({name: joined[name].to_numpy() for name in joined})
    else:
        data = _time_columns([])
    frame = TimeSeriesDataFrame(data)
    frame.attrs = {**provenance, "variables": metadata}
    return frame



def _clip_variable_dict(data, trange, *, layouts=None):
    """Clip decoded backend variables to an inclusive time interval."""
    start, end = normalize_interval(trange)
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


def _merge_variable_dicts(chunks, *, layouts=None):
    """Merge decoded backend records while preserving their native axes."""
    grouped = {}
    for chunk in chunks:
        for name, variable in chunk.items():
            grouped.setdefault(name, []).append(variable)
    merged = {}
    for name, parts in grouped.items():
        dates = pd.DatetimeIndex(
            np.concatenate([Time(part["x"]).datetime.to_numpy() for part in parts])
        )
        order = np.argsort(dates.asi8, kind="stable")
        order = order[~dates[order].duplicated(keep="first")]
        item = deepcopy(parts[0])
        item["x"] = dates[order]
        layout = (layouts or {}).get(name, {})
        time_axis = layout.get("time_axis", 0)
        try:
            item["y"] = np.take(
                np.concatenate([np.asarray(part["y"]) for part in parts], axis=time_axis),
                order,
                axis=time_axis,
            )
        except ValueError as exc:
            raise ValueError(
                f"Variable {name!r} changes sample shape between intervals"
            ) from exc
        axis_keys = {
            key for part in parts for key in part if re.fullmatch(r"v\d*", key)
        }
        for key in axis_keys:
            if any(key not in part for part in parts):
                raise ValueError(
                    f"Variable {name!r} is missing axis {key!r} in an interval"
                )
            axes = [np.asarray(part[key]) for part in parts]
            time_axes = [
                layout.get("axes", {}).get(key, {}).get(
                    "time_axis", None if axis.ndim <= 1 else 0
                )
                for axis in axes
            ]
            if all(axis_time is None for axis_time in time_axes) and all(
                np.array_equal(axis, axes[0]) for axis in axes
            ):
                item[key] = axes[0].copy()
            else:
                samples = [
                    np.broadcast_to(axis, (len(part["x"]), *axis.shape))
                    if axis_time is None
                    else np.moveaxis(axis, axis_time, 0)
                    for part, axis, axis_time in zip(parts, axes, time_axes)
                ]
                combined = np.concatenate(samples, axis=0)[order]
                item[key] = (
                    np.moveaxis(combined, 0, time_axes[0])
                    if time_axes[0] is not None
                    else combined
                )
        merged[name] = item
    return merged

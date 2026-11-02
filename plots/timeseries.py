"""Generic time-series plots and their variable-selection helpers."""

import numpy as np
import pandas as pd

from common.frame import read_variable
from plots.labels import VariableInfo, axis_label, describe
from plots.tools import pass_axes

__all__ = ["tplot"]


def _available_variables(data):
    if not isinstance(data, pd.DataFrame):
        return list(data)
    metadata = data.attrs.get("variables", {})
    available = [
        name for name, info in metadata.items()
        if all(column in data.columns for column in info.get("columns", [name]))
    ]
    if available:
        return available
    return [name for name in data.columns if name not in ("Datetime", "Timestamp")]


def resolve_dataset_key(data, key, prefix=None):
    """Find the full variable name requested by a plotting caller.

    Mission datasets often store qualified names such as
    ``thb_peir_density``, while callers may provide either that full name or the
    shorter ``peir_density``.  Resolution proceeds in this order:

    1. Return ``key`` when it is already an exact variable name.
    2. When ``prefix`` is supplied, try ``<prefix>_<key>``; a trailing
       underscore on ``prefix`` is removed before joining.
    3. Otherwise, accept a unique variable whose name ends in ``_<key>``.

    A suffix shared by multiple variables is ambiguous and raises ``KeyError``
    instead of silently choosing one.  A missing variable also raises
    ``KeyError`` and includes a sample of the available names.  This helper lets
    :func:`tplot` accept convenient short names without losing deterministic
    selection when a dataset contains multiple probes or missions.
    """
    keys = _available_variables(data)
    if key in keys:
        return key
    if prefix is not None:
        prefixed = f"{str(prefix).rstrip('_')}_{key}"
        if prefixed in keys:
            return prefixed
    matches = sorted(name for name in keys if name.endswith("_" + key))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise KeyError(f"Ambiguous variable {key!r}; matches: {matches}")
    raise KeyError(f"Cannot resolve variable {key!r}; available examples: {keys[:8]}")


def dataset_variable(data, key, prefix=None):
    """Return one variable in legacy ``x/y/v`` form from either contract."""
    resolved = resolve_dataset_key(data, key, prefix=prefix)
    if isinstance(data, pd.DataFrame):
        metadata = data.attrs.get("variables", {})
        if resolved in metadata:
            variable = read_variable(data, resolved)
            info = metadata[resolved]
            return resolved, {
                "x": variable["time"],
                "y": variable["values"],
                **variable["axes"],
                "label": info.get("label"),
                "units": info.get("units"),
            }
        dates = pd.to_datetime(data["Datetime"] if "Datetime" in data else data["Timestamp"], unit=None if "Datetime" in data else "s")
        values = data[resolved].to_numpy()
        mask = ~pd.isna(values)
        return resolved, {"x": dates[mask], "y": values[mask]}
    return resolved, data[resolved]


@pass_axes()
def tplot( data, key=None, trange=None, prefix=None, f=None, ax=None, label=None, component_labels=None, **kwargs,):
    """Plot one scalar or three-component time-series variable.

    Display labels and units come from ``plots.labels`` unless the caller or
    the loader metadata supplies them. The theme controls styling; this
    function only draws data that the caller has already loaded.
    """
    resolved = None
    info = VariableInfo()
    if key is None and isinstance(data, dict) and {"x", "y"}.issubset(data):
        variable = data
    else:
        if key is None:
            available = _available_variables(data)
            if len(available) != 1:
                raise ValueError(
                    "key is required when the dataset contains multiple variables; "
                    f"available examples: {available[:8]}"
                )
            key = available[0]
        resolved, variable = dataset_variable(data, key, prefix=prefix)
        info = describe(resolved)
        if variable.get("label") or variable.get("units"):
            info = VariableInfo(
                label=variable.get("label") or info.label,
                unit=variable.get("units") or info.unit,
                scale=info.scale,
            )

    if label is None:
        label = variable.get("label") or info.label or resolved or "value"

    t = pd.to_datetime(variable["x"])
    y = np.asarray(variable["y"])
    if trange is not None:
        t0, t1 = pd.to_datetime(trange[0]), pd.to_datetime(trange[1])
        selected_idx = (t >= t0) & (t <= t1)
        t = t[selected_idx]
        y = y[selected_idx]

    if y.ndim == 1:
        ax.plot(t, y, label=label, **kwargs)
    elif y.ndim == 2 and y.shape[1] == 3:
        if component_labels is None:
            base = resolved or label or ""
            component_labels = [f"{base}_{suffix}" if base else suffix for suffix in ("x", "y", "z")]
        colors = ["r", "g", "b"]
        for i, comp_label in enumerate(component_labels):
            line_kwargs = dict(kwargs)
            line_kwargs.setdefault("color", colors[i])
            ax.plot(t, y[:, i], label=comp_label, **line_kwargs)
    else:
        raise ValueError(f"y should be 1-D or shape (N, 3), got {y.shape}")

    # Leave an explicitly supplied ylabel untouched; otherwise fill it from
    # the variable registry so scalar panels are labeled consistently.
    if not ax.get_ylabel():
        ylabel = axis_label(info, resolved)
        if ylabel:
            ax.set_ylabel(ylabel)
    return f, ax

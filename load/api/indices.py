"""Read Kyoto and GFZ indices, declare scalar layouts, and normalize records.

``read_indices`` returns ``(data, layouts)`` like ``read_spedas``. Source
responses are cached by dataset and UTC month under ``prefix/indices``.
``parse_indices`` converts the records to descriptions accepted by assemble.
"""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from common.config import get_configured_paths
from common.time import Time

__all__ = ["INDICES_FIELDS", "read_indices", "parse_indices"]

_KYOTO_URL = "https://wdc.kugi.kyoto-u.ac.jp/hapi/data"
_GFZ_URL = "https://kp.gfz.de/app/json/"
_DEFAULT_PRODUCTS = ("ae", "al", "au")

INDICES_FIELDS = {
    "ae": {"source": "kyoto", "dataset": "min_ae", "field": "aeValue", "units": "nT"},
    "al": {"source": "kyoto", "dataset": "min_ae", "field": "alValue", "units": "nT"},
    "au": {"source": "kyoto", "dataset": "min_ae", "field": "auValue", "units": "nT"},
    "ao": {"source": "kyoto", "dataset": "min_ae", "field": "aoValue", "units": "nT"},
    "dst": {"source": "kyoto", "dataset": "hour_dst", "field": "dstValue", "units": "nT"},
    "sym_h": {"source": "kyoto", "dataset": "min_asysym", "field": "symHValue", "units": "nT"},
    "sym_d": {"source": "kyoto", "dataset": "min_asysym", "field": "symDValue", "units": "nT"},
    "asy_h": {"source": "kyoto", "dataset": "min_asysym", "field": "asyHValue", "units": "nT"},
    "asy_d": {"source": "kyoto", "dataset": "min_asysym", "field": "asyDValue", "units": "nT"},
    "kp": {"source": "gfz", "dataset": "Kp", "field": "Kp", "units": "1"},
    "ap": {"source": "gfz", "dataset": "ap", "field": "ap", "units": "2 nT"},
    "ap_daily": {"source": "gfz", "dataset": "Ap", "field": "Ap", "units": "2 nT"},
    "sn": {"source": "gfz", "dataset": "SN", "field": "SN", "units": "1"},
    "f107": {"source": "gfz", "dataset": "Fobs", "field": "Fobs", "units": "sfu"},
    "f107adj": {"source": "gfz", "dataset": "Fadj", "field": "Fadj", "units": "sfu"},
}
_ALIASES = {"f107obs": "f107", "f10.7": "f107", "f10.7adj": "f107adj"}


def _products(product):
    names = list(_DEFAULT_PRODUCTS) if product is None else [product] if isinstance(product, str) else list(product)
    if not names:
        raise ValueError("product sequence cannot be empty")
    names = [_ALIASES.get(str(name).lower(), str(name).lower()) for name in names]
    unknown = set(names).difference(INDICES_FIELDS)
    if unknown:
        raise ValueError(f"Unsupported index products {sorted(unknown)}; choose from {list(INDICES_FIELDS)}")
    return list(dict.fromkeys(names))


def _interval(trange):
    dates = Time(trange).datetime
    return dates[0], dates[1]


def _cached_json(url, params, path, download, timeout, validate):
    if path.exists() and download != "refresh":
        payload = json.loads(path.read_text())
        validate(payload)
        return payload
    if download in (False, "local", "offline"):
        raise FileNotFoundError(f"Index cache not found: {path}; allow downloads to fetch it")
    import requests

    response = requests.get(url, params=params, timeout=timeout)
    response.raise_for_status()
    payload = response.json()
    validate(payload)
    # Validate before replacing an existing cache; failed refreshes preserve it.
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, suffix=".tmp", delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(payload, stream)
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return payload


def _kyoto_records(payload, names):
    status = payload.get("status", {})
    if status.get("code") not in (1200, 1201):
        raise ValueError(f"Kyoto HAPI error: {status}")
    parameters = payload["parameters"]
    fields = {item["name"]: index for index, item in enumerate(parameters)}
    rows = payload["data"]
    if any(len(row) != len(parameters) for row in rows):
        raise ValueError("Kyoto HAPI rows do not match parameter metadata")
    times = Time([row[fields["Time"]] for row in rows]).datetime
    records = {}
    for name in names:
        spec = INDICES_FIELDS[name]
        index = fields[spec["field"]]
        meta = parameters[index]
        values = np.asarray([row[index] for row in rows], dtype=float)
        if meta.get("fill") is not None:
            values[values == float(meta["fill"])] = np.nan
        attrs = {key: deepcopy(payload[key]) for key in (
            "title", "cadence", "timeStampLocation", "resourceID", "datasetCitation", "modificationDate",
        ) if key in payload}
        attrs.update(source="kyoto", dataset=spec["dataset"], url=_KYOTO_URL, field=spec["field"])
        if "versionCode" in fields:
            attrs["processing_versions_in_source"] = sorted({row[fields["versionCode"]] for row in rows})
        records[name] = {"x": times, "y": values, "units": meta.get("units") or spec["units"], "attrs": attrs}
    return records


def _gfz_records(payload, names):
    if "error" in payload or "datetime" not in payload:
        raise ValueError(f"GFZ response error: {payload}")
    times = Time(payload["datetime"]).datetime
    records = {}
    for name in names:
        spec = INDICES_FIELDS[name]
        values = np.asarray(payload[spec["field"]], dtype=float)
        if values.ndim != 1 or len(values) != len(times):
            raise ValueError(f"GFZ {spec['field']} values do not match timestamps")
        # GFZ marks unavailable geomagnetic and solar index values as negative.
        values[values < 0] = np.nan
        attrs = deepcopy(payload.get("meta", {}))
        attrs.update(source="gfz", dataset=spec["dataset"], url=_GFZ_URL, field=spec["field"], timeStampLocation="begin")
        if "status" in payload:
            if len(payload["status"]) != len(times):
                raise ValueError("GFZ status does not match timestamps")
            attrs["processing_statuses_in_source"] = sorted(set(payload["status"]))
        records[name] = {"x": times, "y": values, "units": spec["units"], "attrs": attrs}
    return records


def read_indices(trange, product=None, *, prefix=None, download="local_first", timeout=30):
    """Return ``(data, layouts)`` for one interval from Kyoto/GFZ.

    Products default to AE/AL/AU. ``prefix`` is the cache root (otherwise the
    configured SPEDAS directory). False/local/offline uses cache only;
    True/download/local_first/auto fetches missing months; refresh replaces
    cached months. Source timestamps are retained, including Kyoto Dst at
    HH:29:30 and GFZ daily date labels at 00:00 UTC. Missing source values
    become NaN; failed requests and invalid responses raise, never fake data.
    """
    names = _products(product)
    start, end = _interval(trange)
    if download not in (False, True, "local", "offline", "download", "local_first", "auto", "refresh"):
        raise ValueError("download must be local/offline, download/local_first/auto, or refresh")
    root = Path(prefix if prefix is not None else get_configured_paths()["spedas_base"]).expanduser() / "indices"
    groups = {}
    for name in names:
        spec = INDICES_FIELDS[name]
        groups.setdefault((spec["source"], spec["dataset"]), []).append(name)
    chunks = {name: [] for name in names}
    for month in pd.period_range(start.to_period("M"), end.to_period("M"), freq="M"):
        first, stop = month.start_time, (month + 1).start_time
        for (source, dataset), selected in groups.items():
            if source == "kyoto":
                url = _KYOTO_URL
                params = {"id": dataset, "time.min": first.isoformat() + "Z", "time.max": stop.isoformat() + "Z", "format": "json", "include": "header"}
                decode = _kyoto_records
            else:
                url = _GFZ_URL
                params = {"index": dataset, "start": first.isoformat() + "Z", "end": stop.isoformat() + "Z"}
                decode = _gfz_records
            path = root / source / dataset / f"{month.strftime('%Y%m')}.json"
            payload = _cached_json(url, params, path, download, timeout, lambda data: decode(data, selected))
            records = decode(payload, selected)
            for name, item in records.items():
                times = item["x"]
                mask = np.asarray((times >= start) & (times <= end) & (times >= first) & (times < stop))
                chunks[name].append({**item, "x": times[mask], "y": item["y"][mask]})
    data = {}
    for name, parts in chunks.items():
        dates = pd.DatetimeIndex(np.concatenate([part["x"].to_numpy() for part in parts]))
        values = np.concatenate([part["y"] for part in parts])
        order = np.argsort(dates.asi8, kind="stable")
        order = order[~dates[order].duplicated(keep="first")]
        attrs = deepcopy(parts[0]["attrs"])
        # Keep month-specific metadata without silently discarding later versions.
        attrs["source_chunks"] = [deepcopy(part["attrs"]) for part in parts]
        data[name] = {"x": dates[order], "y": values[order], "units": parts[0]["units"], "attrs": attrs}
    layouts = {name: {"kind": "scalar", "time_axis": 0, "columns": [name]} for name in names}
    return data, layouts


def parse_indices(data, layouts):
    """Normalize decoded index records into explicit scalar variable descriptions."""
    variables = {}
    for name, item in data.items():
        layout = layouts[name]
        values = np.asarray(item["y"], dtype=float)
        if layout["kind"] != "scalar" or layout["time_axis"] != 0 or len(layout["columns"]) != 1:
            raise ValueError(f"Index {name!r} requires a scalar layout")
        if values.ndim != 1 or len(values) != len(item["x"]):
            raise ValueError(f"Index {name!r} values do not match timestamps")
        variables[name] = {
            "time": item["x"], "values": values, "time_axis": 0,
            "kind": "scalar", "columns": list(layout["columns"]),
            "units": item.get("units"), "attrs": deepcopy(item.get("attrs", {})),
        }
    return variables

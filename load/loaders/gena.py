
from __future__ import annotations

from pathlib import Path
import re

import pandas as pd
import pyarrow.parquet as pq

from common.config import get_configured_paths
from gena.frame import ParameterDataFrame, ParticleDataFrame
from common.time import Time


#: Columns that only engineering/PDF rows carry.
_ENGINEERING_MARKERS = {"MJ2000_Star_X", "MJ2000_Star_Y", "MJ2000_Star_Z"}
#: Columns that identify particle-event rows, including L1A TDC channels.
_EVENT_MARKERS = {"TX_1s", "TX_2s", "TY_1s", "TY_2s", "X1", "Y1", "Tof", "EP"}


def _frame_type(product, columns):
    """Choose the frame class from the declared product and the actual columns."""
    if product == "pdf":
        return ParameterDataFrame
    present = set(columns)
    if _ENGINEERING_MARKERS <= present and not (_EVENT_MARKERS & present):
        return ParameterDataFrame
    return ParticleDataFrame


def _provenance(product):
    """Split a GENA product suffix into level and probe metadata."""
    if not product:
        return {"source": "gena", "product": None, "level": None, "probe": None}
    if product[-2:] in {"_a", "_b"}:
        return {
            "source": "gena",
            "product": product,
            "level": product[:-2],
            "probe": product[-1],
        }
    return {"source": "gena", "product": product, "level": product, "probe": None}




def load_gena(*, file=None, trange=None, product=None, columns=None):
    """Load one complete GENA Parquet file or a product over a time range."""
    if (file is None) == (trange is None):
        raise ValueError("provide exactly one of file or trange")
    if file is not None and product is not None:
        raise ValueError("product cannot be used with file")
    if trange is not None and product is None:
        raise ValueError("product is required with trange")
    if columns is not None:
        columns = list(dict.fromkeys(["Timestamp", *columns]))
    start = end = None
    if trange is not None:
        start, end = Time(trange).datetime
    inferred_product = None
    if file is not None:
        file_paths = [Path(file).expanduser()]
        if not file_paths[0].is_file():
            raise FileNotFoundError(f"Parquet file does not exist: {file_paths[0]}")
        match = re.fullmatch(r"\d{8}_(.+)\.pqt", file_paths[0].name)
        if match:
            inferred_product = match.group(1)
    else:
        root = Path(get_configured_paths()["pqt_base"])
        filename_pattern = re.compile(rf"(\d{{8}})_{re.escape(product)}\.pqt")
        requested_days = {
            day.strftime("%Y%m%d")
            for day in pd.date_range(start.floor("D"), end.floor("D"), freq="D")
        }
        file_paths = []
        if root.is_dir():
            for candidate in root.rglob(f"*_{product}.pqt"):
                match = filename_pattern.fullmatch(candidate.name)
                if match and match.group(1) in requested_days:
                    file_paths.append(candidate)
        file_paths.sort()
    effective_product = product or inferred_product
    filters = None
    if trange is not None:
        filters = [
            ("Timestamp", ">=", start.timestamp()),
            ("Timestamp", "<=", end.timestamp()),
        ]
    frames = []
    schema_columns = None
    for file_path in file_paths:
        current_columns = pq.ParquetFile(file_path).schema_arrow.names
        if "Timestamp" not in current_columns:
            raise ValueError(f"Timestamp column is missing from {file_path}")
        if columns is not None:
            missing = [name for name in columns if name not in current_columns]
            if missing:
                raise KeyError(f"columns not found in Parquet data: {missing}")
        if schema_columns is None:
            schema_columns = current_columns
        read_columns = columns if effective_product == "pdf" else None
        frames.append(
            pq.read_table(file_path, filters=filters, columns=read_columns).to_pandas()
        )
    result = (
        pd.concat(frames, ignore_index=True)
        if frames
        else pd.DataFrame(columns=list(columns or ["Timestamp"]))
    )
    if trange is not None and not result.empty:
        result = result.loc[
            result["Timestamp"].between(
                start.timestamp(), end.timestamp(), inclusive="both"
            )
        ]
    if not result.empty:
        if effective_product == "pdf":
            result = result.drop_duplicates(subset=["Timestamp"], keep="first")
        else:
            # Event identity is (UUID, PID). PacketID/EventID are per-file
            # provenance indexes and must not participate in global dedupe.
            if {"UUID", "PID"}.issubset(result.columns):
                result = result.drop_duplicates(subset=["UUID", "PID"], keep="first")
            else:
                event_columns = [name for name in result.columns if name != "UUID"]
                result = result.drop_duplicates(subset=event_columns, keep="first")
        result = result.sort_values("Timestamp", kind="stable")
    if columns is not None:
        missing = [name for name in columns if name not in result.columns]
        if missing and schema_columns is not None:
            raise KeyError(f"columns not found in Parquet data: {missing}")
        result = result.reindex(columns=columns)
    result = result.reset_index(drop=True)
    result["Timestamp"] = pd.to_numeric(result["Timestamp"], errors="raise").astype(float)
    result["Datetime"] = pd.to_datetime(
        result["Timestamp"].to_numpy(dtype=float), unit="s", utc=True
    ).tz_localize(None)
    frame = _frame_type(effective_product, result.columns)(result)
    frame.attrs.update(_provenance(effective_product))
    return frame





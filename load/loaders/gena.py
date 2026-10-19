
from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from common.config import get_configured_paths
from gena.frame import ParameterDataFrame, ParticleDataFrame
from load.api.public import normalize_interval




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
    intervals = []
    if trange is not None:
        intervals = [tuple(normalize_interval(trange))]
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
            for start, end in intervals
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
    if intervals:
        filters = [
            [
                ("Timestamp", ">=", start.timestamp()),
                ("Timestamp", "<=", end.timestamp()),
            ]
            for start, end in intervals
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
    if intervals and not result.empty:
        mask = np.zeros(len(result), dtype=bool)
        for start, end in intervals:
            mask |= result["Timestamp"].between(
                start.timestamp(), end.timestamp(), inclusive="both"
            )
        result = result.loc[mask]
    if not result.empty:
        if effective_product == "pdf":
            result = result.drop_duplicates(subset=["Timestamp"], keep="first")
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
    frame_type = ParameterDataFrame if effective_product == "pdf" else ParticleDataFrame
    return frame_type(result)





"""Schema-checked, overwrite-based Parquet output for raw conversion."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from uuid import uuid4
import os

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from . import structure as st

PARQUET_ENGINE = "pyarrow"
COMPRESSION = "zstd"
SCHEMA_VERSION = 1

SUPPORTED_SUFFIX = {
    "pdf": st.pdf_type,
    "l1_a": st.df1_type,
    "l1_b": st.df1_type,
    "l1a_a": st.df1a_type,
    "l1a_b": st.df1a_type,
    "l2_a": st.df2_type,
    "l2_b": st.df2_type,
    "l2l_a": st.df2_type,
    "l2l_b": st.df2_type,
    "l2m_a": st.df2_type,
    "l2m_b": st.df2_type,
    "l2h_a": st.df2_type,
    "l2h_b": st.df2_type,
    "l2f_a": st.df2_type,
    "l2f_b": st.df2_type,
}


def _pa_type(dtype):
    dtype = np.dtype(dtype)
    if dtype.kind == "S":
        return pa.string()
    if dtype.kind == "b":
        return pa.bool_()
    if dtype.kind == "u":
        return pa.uint8() if dtype.itemsize == 1 else (
            pa.uint16() if dtype.itemsize == 2 else (
                pa.uint32() if dtype.itemsize == 4 else pa.uint64()))
    if dtype.kind == "i":
        return pa.int8() if dtype.itemsize == 1 else (
            pa.int16() if dtype.itemsize == 2 else (
                pa.int32() if dtype.itemsize == 4 else pa.int64()))
    if dtype.kind == "f":
        return pa.float32() if dtype.itemsize == 4 else pa.float64()
    raise TypeError(f"unsupported Parquet dtype: {dtype}")


def product_schema(suffix, metadata=None):
    """Return the Arrow schema for one product suffix."""
    if suffix not in SUPPORTED_SUFFIX:
        raise ValueError(f"unsupported product suffix: {suffix}")
    fields = [(name, _pa_type(dtype)) for name, dtype in SUPPORTED_SUFFIX[suffix].items()]
    return pa.schema(fields, metadata=metadata)


def _cast_series(series, dtype, name):
    dtype = np.dtype(dtype)
    if dtype.kind == "S":
        return series.astype(str)
    if dtype.kind == "b":
        return series.fillna(False).astype(bool)
    if dtype.kind == "f":
        return series.astype(dtype)
    if dtype.kind in "iu":
        if series.isna().any():
            raise ValueError(f"{name} contains NaN but declares {dtype}")
        return series.astype(dtype)
    raise TypeError(f"unsupported dtype for {name}: {dtype}")


def cast_frame(df, suffix):
    """Check columns and cast a frame to the declared product schema."""
    if suffix not in SUPPORTED_SUFFIX:
        raise ValueError(f"unsupported product suffix: {suffix}")
    expected = SUPPORTED_SUFFIX[suffix]
    missing = [name for name in expected if name not in df.columns]
    extra = [name for name in df.columns if name not in expected]
    if missing:
        raise ValueError(f"{suffix} missing columns: {missing}")
    if extra:
        raise ValueError(f"{suffix} unexpected columns: {extra}")
    return pd.DataFrame({name: _cast_series(df[name], dtype, name)
                         for name, dtype in expected.items()})


def _read_metadata(path):
    try:
        metadata = pq.ParquetFile(path).metadata.metadata or {}
    except Exception:
        return {}
    return {key.decode() if isinstance(key, bytes) else str(key):
            value.decode() if isinstance(value, bytes) else str(value)
            for key, value in metadata.items()}


class ProductWriter:
    """Write per-day products to temporary files, then replace atomically.

    ``overwrite=True`` replaces an existing non-empty product. If the existing
    file records a different ``snid`` or schema version, ``force`` is required.
    Temporary files are written under the destination directory and removed on
    abort, so a failed conversion never publishes a partial Parquet file.
    """

    def __init__(self, root, *, by_day=True, overwrite=True, force=False,
                 metadata=None):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.by_day = bool(by_day)
        self.overwrite = bool(overwrite)
        self.force = bool(force)
        self.metadata = {str(k).encode(): str(v).encode()
                         for k, v in (metadata or {}).items()}
        self._writers = {}
        self._temps = {}
        self._schemas = {}
        self._rows = Counter()
        self._paths = {}

    def _final_path(self, day, suffix):
        name = f"{day}_{suffix}.pqt" if day else f"{suffix}.pqt"
        return self.root / name

    def _check_existing(self, final):
        if not final.exists():
            return
        if not self.overwrite:
            raise FileExistsError(f"refusing to overwrite {final}")
        metadata = _read_metadata(final)
        existing_snid = metadata.get("gena_snid")
        if existing_snid is not None and str(existing_snid) != str(self.metadata.get(b"gena_snid", b"").decode()):
            if not self.force:
                raise ValueError(
                    f"{final} was written with snid={existing_snid}; "
                    "use force=True to overwrite a different calibration")
        existing_version = metadata.get("gena_schema_version")
        if existing_version is not None and str(existing_version) != str(SCHEMA_VERSION):
            if not self.force:
                raise ValueError(
                    f"{final} uses schema version {existing_version}; "
                    "use force=True to overwrite")

    def _open(self, day, suffix):
        key = (day, suffix)
        if key in self._writers:
            return self._writers[key], self._schemas[key]
        final = self._final_path(day, suffix)
        self._check_existing(final)
        temp = final.with_name(f".{final.name}.{os.getpid()}.{uuid4().hex}.tmp")
        schema = product_schema(suffix, metadata=self.metadata)
        writer = pq.ParquetWriter(temp, schema, compression=COMPRESSION)
        self._writers[key] = writer
        self._schemas[key] = schema
        self._temps[key] = temp
        self._paths[key] = final
        return writer, schema

    def write(self, df, suffix):
        """Append a frame to the temporary product for its UTC day."""
        if df.empty:
            return
        frame = cast_frame(df, suffix)
        if self.by_day:
            if "Timestamp" not in frame.columns:
                raise ValueError(f"{suffix} requires Timestamp for daily output")
            dates = pd.to_datetime(frame["Timestamp"], unit="s", utc=True).dt.strftime("%Y%m%d")
            if dates.isna().any():
                raise ValueError(f"{suffix} contains invalid Timestamp values")
            groups = frame.groupby(dates, sort=True)
        else:
            groups = [("", frame)]
        for day, part in groups:
            writer, schema = self._open(day, suffix)
            table = pa.Table.from_pandas(
                part.reset_index(drop=True), schema=schema,
                preserve_index=False, safe=True)
            writer.write_table(table)
            self._rows[(day, suffix)] += len(part)

    def finish(self):
        """Close all writers and atomically replace the final products."""
        try:
            for writer in self._writers.values():
                writer.close()
            self._writers.clear()
            for key, temp in list(self._temps.items()):
                final = self._paths[key]
                if self._rows[key] == 0:
                    temp.unlink(missing_ok=True)
                    continue
                self._check_existing(final)
            for key, temp in list(self._temps.items()):
                if self._rows[key] == 0:
                    continue
                os.replace(temp, self._paths[key])
            rows = {self._paths[key].name: int(self._rows[key]) for key in self._paths}
            paths = {self._paths[key].name: self._paths[key] for key in self._paths}
            return rows, paths
        except BaseException:
            self.abort()
            raise
        finally:
            self._temps.clear()
            self._paths.clear()

    def abort(self):
        for writer in self._writers.values():
            try:
                writer.close()
            except Exception:
                pass
        self._writers.clear()
        for temp in self._temps.values():
            try:
                temp.unlink(missing_ok=True)
            except OSError:
                pass
        self._temps.clear()
        self._paths.clear()

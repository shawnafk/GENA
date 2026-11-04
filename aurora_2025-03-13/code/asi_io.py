"""Minimal reader for the THEMIS ASI PGM multi-image files used by the UCalgary archive.

A file holds one minute of 256x256 16-bit frames (3 s cadence).  Every frame is a
complete P5 (binary PGM) image with its own comment header that carries the site,
the site coordinates, the optical projection and the exposure start time.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


@dataclass
class Frame:
    time: datetime
    image: np.ndarray
    meta: dict


def _read_lines(buf: bytes, pos: int):
    end = buf.index(b"\n", pos)
    return buf[pos:end].decode("ascii", "replace"), end + 1


def _header_dict(lines: list[str]) -> dict:
    meta: dict = {}
    for line in lines:
        if not line.startswith("#"):
            continue
        body = line.lstrip("#").strip()
        if body.startswith('"'):
            parts = body.split('"')
            if len(parts) >= 3:
                meta[parts[1]] = parts[2].strip()
        elif " " in body:
            key, val = body.split(" ", 1)
            meta[key] = val.strip()
    return meta


def read_frames(path: str | Path, max_frames: int | None = None) -> list[Frame]:
    path = Path(path)
    raw = gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
    frames: list[Frame] = []
    pos = 0
    n = len(raw)
    while pos < n:
        if raw[pos : pos + 3] != b"P5\n":
            raise ValueError(f"unexpected magic at byte {pos}: {raw[pos:pos+8]!r}")
        pos += 3
        comments: list[str] = []
        while raw[pos : pos + 1] == b"#":
            line, pos = _read_lines(raw, pos)
            comments.append(line)
        dims_line, pos = _read_lines(raw, pos)
        width, height = (int(v) for v in dims_line.split())
        maxval_line, pos = _read_lines(raw, pos)
        maxval = int(maxval_line)
        nbytes = width * height * (1 if maxval < 256 else 2)
        chunk = raw[pos : pos + nbytes]
        pos += nbytes
        dtype = ">u2" if maxval >= 256 else np.uint8
        image = np.frombuffer(chunk, dtype=dtype).reshape(height, width).astype(np.float32)
        meta = _header_dict(comments)
        stamp = meta.get("Image request start", "")
        when = None
        if stamp:
            head, _, frac = stamp.partition(".")
            fmt = "%Y-%m-%d %H:%M:%S"
            base = datetime.strptime(head.strip(), fmt)
            micro = int((frac.split()[0] + "000000")[:6]) if frac else 0
            when = base.replace(microsecond=micro, tzinfo=timezone.utc)
        frames.append(Frame(when, image, meta))
        if max_frames is not None and len(frames) >= max_frames:
            break
    return frames


def read_row_series(path: str | Path) -> tuple[np.ndarray, list[datetime], dict]:
    """Reduce a `row2` file (752/2 x 1 pixel strip every 3 s) to a (time, x) image."""
    frames = read_frames(path)
    rows = np.vstack([f.image.reshape(-1) for f in frames])
    times = [f.time for f in frames]
    return rows, times, frames[0].meta


def circle_mask(size: int = 256, center: tuple[float, float] | None = None, radius: float | None = None):
    cy, cx = center if center else ((size - 1) / 2, (size - 1) / 2)
    r = radius if radius is not None else size / 2
    yy, xx = np.mgrid[0:size, 0:size]
    return (yy - cy) ** 2 + (xx - cx) ** 2 <= r**2

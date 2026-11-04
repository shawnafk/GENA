"""Render THEMIS ASI frames, per-station sequences, a station overview and keograms.

Display recipe (same for every panel so frames are comparable):
  * 3x3 median filter removes stars and cosmic-ray hits,
  * sky level = median of the circular field (r = 0.97 R),
  * linear stretch to the 99.7th percentile of the field, gamma 0.7,
  * everything outside the all-sky circle is black.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402
from scipy.ndimage import median_filter  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from asi_io import circle_mask, read_frames  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "_raw"
SITE_ORDER = ["atha", "fsim", "fsmi", "gako", "gill", "inuv",
              "kapu", "pina", "rank", "snkq", "talo", "tpas"]
MASK = circle_mask(256, radius=124)
OUTSIDE = ~MASK


def load_display(path: Path, frame_index: int = -1):
    frames = read_frames(path)
    idx = frame_index if frame_index >= 0 else len(frames) // 2
    frame = frames[min(idx, len(frames) - 1)]
    site_meta = frames[0].meta  # only the first frame of a file carries the site block
    img = median_filter(frame.image, 3)
    sky = float(np.median(img[MASK]))
    data = img - sky
    vmax = float(np.percentile(data[MASK], 99.7))
    vmax = max(vmax, 1.0)
    scaled = np.clip(data / vmax, 0.0, 1.0) ** 0.7
    scaled[OUTSIDE] = 0.0
    return scaled, frame, site_meta


def site_label(meta: dict) -> str:
    name = meta.get("Site name", "").strip()
    lat = float(meta.get("Geodetic latitude", "nan"))
    lon = float(meta.get("Geodetic longitude", meta.get("Geodetic Longitude", "nan")))
    return f"{name} ({abs(lat):.1f}\u00b0{'N' if lat >= 0 else 'S'}, {abs(lon):.1f}\u00b0{'E' if lon >= 0 else 'W'})"


def frame_png(path: Path, dest: Path, title_site: str) -> None:
    scaled, frame, meta = load_display(path)
    fig, ax = plt.subplots(figsize=(4.6, 4.95), dpi=170)
    ax.imshow(scaled, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
    ax.add_patch(Circle((127.5, 127.5), 127.5, fill=False, lw=0.6, ec="0.55"))
    ax.set_axis_off()
    ax.set_title(f"{title_site.upper()}  \u2022  {site_label(meta)}", fontsize=7.5, pad=3)
    fig.text(0.5, 0.035, f"{frame.time:%Y-%m-%d %H:%M:%S} UT   \u2022   THEMIS ASI, 3 s exposure",
             ha="center", fontsize=7, color="0.25")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.06)
    fig.savefig(dest, facecolor="white")
    plt.close(fig)


def sequence_png(paths: list[Path], site: str, dest: Path, ncols: int = 4) -> None:
    nrows = int(np.ceil(len(paths) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.1 * ncols, 3.25 * nrows), dpi=130)
    axes = np.atleast_1d(axes).ravel()
    meta = None
    for ax, path in zip(axes, paths):
        scaled, frame, meta = load_display(path)
        ax.imshow(scaled, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
        ax.set_title(f"{frame.time:%H:%M} UT", fontsize=8, pad=2)
        ax.set_axis_off()
    for ax in axes[len(paths):]:
        ax.set_axis_off()
    label = site_label(meta) if meta else site
    fig.suptitle(f"THEMIS ASI {site.upper()}  \u2022  {label}  \u2022  2025-03-13 (UTC night)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(dest, facecolor="white")
    plt.close(fig)


def overview_png(paths: list[Path], dest: Path, stamp: str) -> None:
    ncols, nrows = 4, int(np.ceil(len(paths) / 4))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.15 * ncols, 3.35 * nrows), dpi=140)
    axes = np.atleast_1d(axes).ravel()
    for ax, path in zip(axes, paths):
        scaled, frame, meta = load_display(path)
        ax.imshow(scaled, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
        sid = meta.get("Site unique ID", "?").strip()
        lat = float(meta.get("Geodetic latitude", "nan"))
        ax.set_title(f"{sid.upper()}  {lat:.1f}\u00b0N", fontsize=9, pad=2)
        ax.set_axis_off()
    for ax in axes[len(paths):]:
        ax.set_axis_off()
    fig.suptitle(f"THEMIS all-sky imagers \u2014 2025-03-13 {stamp} UT", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(dest, facecolor="white")
    plt.close(fig)


def keogram_png(site: str, dest: Path, hours: list[str]) -> None:
    from PIL import Image

    strips, labels = [], []
    for hh in hours:
        hits = sorted((ROOT / "quicklook").glob(f"20250313_{hh}_{site}_*_keogram.jpg"))
        if hits:
            strips.append(np.asarray(Image.open(hits[0]).convert("L")))
            labels.append(hh)
    if not strips:
        return
    hour_width = strips[0].shape[1]
    strip = np.hstack(strips)
    fig, ax = plt.subplots(figsize=(13.5, 3.3), dpi=140)
    ax.imshow(strip, cmap="gray", aspect="auto", vmin=0, vmax=255)
    ax.set_xlabel("UT on 2025-03-13", fontsize=8)
    ax.set_ylabel("sky row (keogram)", fontsize=7)
    ax.set_xticks([(i + 0.5) * hour_width for i in range(len(labels))])
    ax.set_xticklabels([f"{h}:00" for h in labels], fontsize=7)
    ax.tick_params(axis="y", labelleft=False, left=False)
    ax.set_title(f"THEMIS ASI {site.upper()} keogram, 2025-03-13", fontsize=9)
    fig.tight_layout()
    fig.savefig(dest, facecolor="white")
    plt.close(fig)


def keogram_stack_png(sites: list[str], dest: Path, hours: list[str]) -> None:
    from PIL import Image

    rows, used = [], []
    for site in sites:
        strips = []
        for hh in hours:
            hits = sorted((ROOT / "quicklook").glob(f"20250313_{hh}_{site}_*_keogram.jpg"))
            strips.append(np.asarray(Image.open(hits[0]).convert("L")) if hits
                          else np.zeros((256, 600), dtype=np.uint8))
        rows.append(np.hstack(strips))
        used.append(site)
    fig, axes = plt.subplots(len(rows), 1, figsize=(14, 1.55 * len(rows)), dpi=150)
    axes = np.atleast_1d(axes)
    for ax, row, site in zip(axes, rows, used):
        ax.imshow(row, aspect="auto", cmap="gray", vmin=0, vmax=255)
        ax.set_ylabel(site.upper(), fontsize=9)
        ax.set_yticks([])
        ax.set_xticks([(i + 0.5) * 600 for i in range(len(hours))])
        ax.set_xticklabels([f"{h}" for h in hours], fontsize=7)
        for i in range(1, len(hours)):
            ax.axvline(i * 600, color="0.6", lw=0.4)
        if site != used[-1]:
            ax.set_xticklabels([])
    axes[0].xaxis.set_ticks_position("top")
    axes[0].tick_params(axis="x", labeltop=True, labelbottom=False, top=True, bottom=False)
    fig.suptitle("THEMIS ASI keograms, 2025-03-13 \u2014 hour labels in UT (blank = no archive file)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(dest, facecolor="white")
    plt.close(fig)


def main() -> int:
    (ROOT / "frames").mkdir(exist_ok=True)
    (ROOT / "mosaics").mkdir(exist_ok=True)
    (ROOT / "overview").mkdir(exist_ok=True)
    (ROOT / "keograms").mkdir(exist_ok=True)

    hours = [f"{h:02d}" for h in range(3, 14)]
    for site in SITE_ORDER:
        paths = sorted(RAW.glob(f"20250313_????_{site}_*_full.pgm.gz"))
        if not paths:
            continue
        outdir = ROOT / "frames" / site
        outdir.mkdir(parents=True, exist_ok=True)
        for path in paths:
            scaled, frame, meta = load_display(path)
            stamp = f"{frame.time:%Y%m%d_%H%M%S}"
            frame_png(path, outdir / f"{site}_{stamp}ut.png", site)
            print("frame", site, stamp, flush=True)
        sequence_png(paths, site, ROOT / "mosaics" / f"{site}_2025-03-13_sequence.png")
        keogram_png(site, ROOT / "keograms" / f"{site}_2025-03-13_keogram.png", hours)
        print("done", site, flush=True)

    overview = sorted(RAW.glob("20250313_0800_*_full.pgm.gz"))
    if overview:
        overview_png(overview, ROOT / "overview" / "stations_2025-03-13_0800ut.png", "08:00")
    return 0


if __name__ == "__main__":
    sys.exit(main())

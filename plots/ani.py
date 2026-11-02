"""GENA time-dependent Matplotlib animations.

The module separates three concerns:

* :func:`gena_fov_frames` turns event data into a reusable
  :class:`FovFrames` snapshot (no Matplotlib);
* :func:`render_fov` encodes those frames to MP4 with progress/cancellation;
* :func:`dump_fov` is the one-call preset used by reports.

``AnimationDemo3`` remains available for callers that constructed it directly;
its original ``run``/``save`` behaviour is preserved.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

from common.config import resolve_output_path
from gena.post import statistic as sts

__all__ = [
    "AnimationDemo3",
    "AnimationSpec",
    "FovFrames",
    "dumpFov",
    "dump_fov",
    "gena_fov_frames",
    "render_fov",
]

ProgressCallback = Callable[[int, int], None]
CancelCallback = Callable[[], bool]


@dataclass(frozen=True)
class AnimationSpec:
    """Encoding options for one animation."""

    fps: int = 5
    dpi: int = 120
    figsize: tuple[float, float] = (10.0, 8.0)
    poster: bool = True
    codec: str = "h264"
    extra_args: tuple[str, ...] = ("-pix_fmt", "yuv420p", "-movflags", "+faststart")

    def __post_init__(self):
        if isinstance(self.fps, bool) or not isinstance(self.fps, (int, np.integer)) or self.fps < 1:
            raise ValueError("fps must be a positive integer")
        if isinstance(self.dpi, bool) or not isinstance(self.dpi, (int, np.integer)) or self.dpi < 1:
            raise ValueError("dpi must be a positive integer")
        if len(self.figsize) != 2 or any(float(value) <= 0 for value in self.figsize):
            raise ValueError("figsize must contain two positive values")
        if not isinstance(self.codec, str) or not self.codec:
            raise ValueError("codec must be a non-empty string")
        if not isinstance(self.extra_args, tuple):
            raise TypeError("extra_args must be a tuple of strings")


@dataclass(frozen=True)
class FovFrames:
    """Precomputed frames for the standard counts/energy/FOV animation."""

    times: np.ndarray
    counts: np.ndarray
    energy_centers: np.ndarray
    energy_counts: np.ndarray
    fov_extent: tuple[float, float, float, float]
    fov_counts: np.ndarray
    cmax: float
    nummax: float

    def __len__(self) -> int:
        return int(len(self.times))


def _frame_indices(count: int, *, max_frames: int | None, stride: int) -> np.ndarray:
    if count < 1:
        return np.empty(0, dtype=int)
    if stride < 1:
        raise ValueError("stride must be a positive integer")
    if max_frames is not None:
        if isinstance(max_frames, bool) or not isinstance(max_frames, (int, np.integer)):
            raise TypeError("max_frames must be a positive integer or None")
        if max_frames < 1:
            raise ValueError("max_frames must be a positive integer or None")
    indices = np.arange(0, count, stride, dtype=int)
    if max_frames is not None and len(indices) > max_frames:
        positions = np.linspace(0, len(indices) - 1, max_frames).round().astype(int)
        indices = indices[np.unique(positions)]
    return indices


def gena_fov_frames(
    l2,
    *,
    window: str = "5min",
    max_frames: int | None = None,
    stride: int = 1,
) -> FovFrames | None:
    """Precompute the standard GENA counts/energy/FOV frames.

    Histograms are evaluated once per time window, not once per rendered frame.
    ``max_frames`` keeps an evenly spaced subset; ``stride`` takes every n-th
    window before that cap. Returns ``None`` when no usable frames exist.
    """
    time, counts, _ = sts.window_hist0d(l2, window=window)
    energy_time, energy_centers, energy_counts, _ = sts.window_hist1d(
        l2,
        "Energy",
        np.arange(0, 101),
        window=window,
    )
    fov_time, _, _, fov_counts, _ = sts.window_hist2d(
        l2,
        ["PX", "PY"],
        [np.arange(-35, 35, 2.5), np.arange(-20, 20, 2.5)],
        window=window,
    )
    frame_count = min(len(time), len(energy_time), len(fov_time))
    if frame_count == 0:
        return None

    indices = _frame_indices(frame_count, max_frames=max_frames, stride=stride)
    if len(indices) == 0:
        return None

    selected_energy = energy_counts[indices]
    selected_fov = fov_counts[indices]
    cmax = max(1.0, float(np.nanpercentile(selected_fov, 99)))
    nummax = max(1.0, float(np.nanpercentile(selected_energy, 99)))
    return FovFrames(
        times=np.asarray(time[:frame_count])[indices],
        counts=np.asarray(counts[:frame_count])[indices],
        energy_centers=np.asarray(energy_centers),
        energy_counts=selected_energy,
        fov_extent=(-35.0, 35.0, -20.0, 20.0),
        fov_counts=selected_fov,
        cmax=cmax,
        nummax=nummax,
    )


def _ffmpeg_path() -> str | None:
    """Return an ffmpeg executable, using imageio-ffmpeg when needed."""
    if animation.writers.is_available("ffmpeg"):
        return None
    from imageio_ffmpeg import get_ffmpeg_exe

    return get_ffmpeg_exe()


def render_fov(
    frames: FovFrames | None,
    output_path,
    *,
    spec: AnimationSpec | None = None,
    progress: ProgressCallback | None = None,
    cancel: CancelCallback | None = None,
) -> Path | None:
    """Encode precomputed FOV frames to MP4.

    Writes to a temporary file and atomically replaces the target, so an
    interrupted run never leaves a truncated MP4. Returns ``None`` for empty
    frames or a cancelled render. A poster JPG is written beside the MP4 when
    ``spec.poster`` is true.
    """
    if frames is None or len(frames) == 0:
        return None
    spec = spec or AnimationSpec()
    destination = resolve_output_path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.stem}.tmp.mp4")
    poster = destination.with_suffix(".jpg")
    demo = AnimationDemo3(
        [frames.times, frames.counts],
        frames.times,
        frames.energy_counts,
        frames.energy_centers,
        frames.fov_counts,
        frames.fov_extent,
        cmax=frames.cmax,
        nummax=frames.nummax,
        figsize=spec.figsize,
    )
    writer = animation.FFMpegWriter(
        fps=spec.fps,
        codec=spec.codec,
        extra_args=list(spec.extra_args),
    )
    completed = False
    try:
        if spec.poster:
            demo.update(0)
            demo.fig.savefig(poster, dpi=spec.dpi, bbox_inches="tight")

        ffmpeg = _ffmpeg_path()
        options = {"animation.ffmpeg_path": ffmpeg} if ffmpeg else {}
        with mpl.rc_context(options):
            with writer.saving(demo.fig, str(temporary), dpi=spec.dpi):
                for index in range(len(frames)):
                    if cancel is not None and cancel():
                        return None
                    demo.update(index)
                    writer.grab_frame()
                    if progress is not None:
                        progress(index + 1, len(frames))
        temporary.replace(destination)
        completed = True
        return destination
    finally:
        plt.close(demo.fig)
        if not completed:
            temporary.unlink(missing_ok=True)


def dump_fov(
    l2,
    output_path,
    *,
    window: str = "5min",
    fps: int = 5,
    spec: AnimationSpec | None = None,
    max_frames: int | None = None,
    progress: ProgressCallback | None = None,
    cancel: CancelCallback | None = None,
) -> Path | None:
    """Build and encode the standard GENA counts/energy/FOV animation."""
    frames = gena_fov_frames(l2, window=window, max_frames=max_frames)
    if frames is None:
        return None
    spec = spec or AnimationSpec(fps=fps)
    return render_fov(frames, output_path, spec=spec, progress=progress, cancel=cancel)


# Compatibility name used by ``gena.post.export`` and older scripts.
dumpFov = dump_fov


class AnimationDemo3:
    """Legacy three-panel animation renderer.

    New code should use :func:`gena_fov_frames` + :func:`render_fov`. This class
    keeps the original constructor and ``run``/``save`` methods available.
    """

    def __init__(
        self,
        tflux,
        t_window,
        energy,
        erange,
        fov_xy,
        xyextent,
        cmax=10,
        nummax=200,
        figsize=(10, 8),
    ):
        self.tflux = tflux[0]
        self.flux = tflux[1]
        self.t_window = t_window
        self.erange = erange
        self.ary1d = energy
        self.ary2d = fov_xy
        self.xyextent = xyextent
        self.nummax = nummax
        self.cmax = cmax

        self.fig, (self.ax1, self.ax2, self.ax3) = plt.subplots(
            3,
            1,
            gridspec_kw={"height_ratios": [1, 1, 3]},
            figsize=figsize,
        )
        self.fig.subplots_adjust(hspace=0.4)

        self.ax1.plot(self.tflux, self.flux, color="blue")
        self.ax1.set_xlabel("Time")
        self.ax1.set_ylabel("Counts")
        self.ax1.grid(True)
        self.ax1.set_title("ENA counts")

        (self.line2,) = self.ax2.plot(
            self.erange,
            self.ary1d[0, :],
            color="blue",
        )
        self.ax2.set_ylabel("Number")
        self.ax2.set_xlabel("E (keV)")
        self.ax2.set_xscale("log")
        self.ax2.set_ylim(0, self.nummax)

        self.frame_id = self.ax1.text(
            0.1,
            0.7,
            "frame = 0",
            transform=self.ax1.transAxes,
        )
        self.im = self.ax3.imshow(
            self.ary2d[0, :, :] * 0,
            aspect="auto",
            extent=self.xyextent,
            vmax=self.cmax,
        )
        self.fig.colorbar(
            self.im,
            ax=self.ax3,
            location="top",
            orientation="horizontal",
            pad=0.1,
        )
        self.current_time_vline = self.ax1.axvline(
            x=self.t_window[0],
            color="r",
            linestyle="--",
        )
        self.ax3.set_aspect("equal")
        self.ax3.set_xlabel("PX (deg)")
        self.ax3.set_ylabel("PY (deg)")
        self.ax3.set_title("Incident-angle distribution")

    def init(self):
        return self.current_time_vline, self.line2, self.im, self.frame_id

    def update(self, frame):
        self.line2.set_data(self.erange, self.ary1d[frame, :])
        self.current_time_vline.set_xdata([self.t_window[frame]])
        self.im.set_array(self.ary2d[frame, :, :])
        self.frame_id.set_text(f"frame = {frame}")
        return self.current_time_vline, self.line2, self.im, self.frame_id

    def run(self, fps=30):
        self.ani = animation.FuncAnimation(
            self.fig,
            self.update,
            frames=len(self.t_window),
            init_func=self.init,
            blit=False,
            interval=1 / fps * 1000,
            repeat=False,
        )
        return self.ani

    def save(self, name, fps=10):
        destination = resolve_output_path(name)
        destination.parent.mkdir(parents=True, exist_ok=True)
        options = {}
        ffmpeg = _ffmpeg_path()
        if ffmpeg:
            options["animation.ffmpeg_path"] = ffmpeg
        with mpl.rc_context(options):
            self.ani.save(destination, writer=animation.FFMpegWriter(fps=fps))

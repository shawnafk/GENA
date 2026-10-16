"""GENA time-dependent Matplotlib animations."""

from __future__ import annotations

from pathlib import Path
from common.config import resolve_output_path

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib import animation

from gena.post import statistic as sts

__all__ = ["AnimationDemo3", "dumpFov"]

## need to fix:
#1. the scale vmax and vmin, for some case it is too small, and the colorbar is not good.
#2. the energy bin should be set to log scale, and the counts range should be in log scale as well, in linear case the line is overflowed since we use small range. 
class AnimationDemo3:
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
            figsize=(10, 8),
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
        if not animation.writers.is_available("ffmpeg"):
            from imageio_ffmpeg import get_ffmpeg_exe

            options["animation.ffmpeg_path"] = get_ffmpeg_exe()
        with mpl.rc_context(options):
            self.ani.save(destination, writer=animation.FFMpegWriter(fps=fps))


def dumpFov(l2, output_path, *, window="5min", fps=5):
    """Create and save the standard GENA count/energy/FOV animation."""
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
    if len(time) == 0 or len(energy_time) == 0 or len(fov_time) == 0:
        return None

    frame_count = min(len(time), len(energy_time), len(fov_time))
    energy_counts = energy_counts[:frame_count]
    fov_counts = fov_counts[:frame_count]
    demo = AnimationDemo3(
        [time[:frame_count], counts[:frame_count]],
        energy_time[:frame_count],
        energy_counts,
        energy_centers,
        fov_counts,
        [-35, 35, -20, 20],
        cmax=max(1, float(np.nanpercentile(fov_counts, 99))),
        nummax=max(1, float(np.nanpercentile(energy_counts, 99))),
    )
    destination = resolve_output_path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        demo.run(fps=fps)
        demo.save(destination, fps=fps)
    finally:
        plt.close(demo.fig)
    return destination




## todo
'''
need to add another animation function for analysis the projection of other object in the sky that enters the fov, such the moon, themis, mms and the Earth itself and Bowshock and Magnetopause, and the plasma sheet. 
'''

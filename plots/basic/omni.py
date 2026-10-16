"""OMNI and geomagnetic-index product panels."""

import matplotlib.pyplot as plt

from plots import vis
from .timeseries import tplot as _tplot

plt.rcParams["legend.frameon"] = False
plt.rcParams["font.size"] = 14


PANEL_CONFIGS = {
    "AE_INDEX": {"ylabel": r"$AE~(nT)$", "color": "k", "zero": True},
    "SYM_H": {"ylabel": r"$SYM_H~(nT)$", "color": "k", "zero": True},
    "BZ_GSM": {"ylabel": r"$B_z~(nT)$", "color": "k", "zero": True},
    "Pressure": {"ylabel": r"$P~(nPa)$", "color": "k"},
    "Beta": {"ylabel": r"$\beta$", "color": "k"},
    "Mgs_mach_num": {"ylabel": r"$M_{mag}$", "color": "k"},
    "Vx": {"ylabel": r"$V_x~(km/s)$", "color": "k"},
    "proton_density": {"ylabel": r"$n_p~(cm^{-3})$", "color": "k"},
}


def plot_panel(data, key, ax, *, trange=None, config=None, **kwargs):
    """Plot one OMNI variable panel from ``data[key]``."""
    cfg = dict(PANEL_CONFIGS.get(key, {}))
    if config:
        cfg.update(config)
    color = kwargs.pop("color", cfg.get("color", None))
    if color is not None:
        kwargs.setdefault("color", color)
    _tplot(data, key=key, trange=trange, f=ax.get_figure(), ax=ax, **kwargs)
    if cfg.get("zero", False):
        ax.axhline(0, ls="--", color="k")
    if "ylabel" in cfg:
        ax.set_ylabel(cfg["ylabel"])
    return ax


def _panel_wrapper(key):
    def wrapper(omni_data, f=None, ax=None, trange=None):
        return plot_panel(omni_data, key, ax, trange=trange)
    return wrapper


_AE = _panel_wrapper("AE_INDEX")
_SYM_H = _panel_wrapper("SYM_H")
_BZ = _panel_wrapper("BZ_GSM")
_Pressure = _panel_wrapper("Pressure")
_Beta = _panel_wrapper("Beta")
_MMN = _panel_wrapper("Mgs_mach_num")
_Vx = _panel_wrapper("Vx")
_Np = _panel_wrapper("proton_density")


@vis.pass_axes(subplot_spec=(6, 1, True, False))
def omni(omni_data, trange=None, keys=None, f=None, ax=None):
    """Plot a standard multi-panel OMNI overview."""
    if keys is None:
        keys = ["BZ_GSM", "Pressure", "Beta", "Mgs_mach_num", "Vx", "proton_density"]
    if len(keys) != len(ax):
        plt.close(f)
        f, axes = plt.subplots(len(keys), 1, sharex=True)
        ax = axes.flatten().tolist()
    for axis, key in zip(ax, keys):
        plot_panel(omni_data, key, axis, trange=trange)
    ax[-1].set_xlabel("Time (UTC)")
    f.autofmt_xdate()
    return f, ax


@vis.pass_axes()
def dst(dst_data, trange=None, f=None, ax=None):
    _tplot(dst_data, trange=trange, f=f, ax=ax, color="k")
    ax.axhline(0, ls="--", color="k")
    return f, ax

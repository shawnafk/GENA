"""OMNI and geomagnetic-index panels built from declarative PanelSpecs."""

from plots import tools
from plots.panels import PanelSpec, plot_panel as _plot_panel, plot_panels

__all__ = ["dst", "omni", "plot_panel"]


OMNI_PANELS = {
    "AE_INDEX": PanelSpec(("AE_INDEX",), zero=True),
    "SYM_H": PanelSpec(("SYM_H",), zero=True),
    "BZ_GSM": PanelSpec(("BZ_GSM",), zero=True),
    "Pressure": PanelSpec(("Pressure",)),
    "Beta": PanelSpec(("Beta",)),
    "Mgs_mach_num": PanelSpec(("Mgs_mach_num",)),
    "Vx": PanelSpec(("Vx",)),
    "proton_density": PanelSpec(("proton_density",)),
}


def plot_panel(data, key, *, trange=None, f=None, ax=None):
    """Plot one OMNI/index variable using its registered PanelSpec."""
    spec = OMNI_PANELS.get(key, PanelSpec((key,)))
    return _plot_panel(data, spec, trange=trange, f=f, ax=ax)


def omni(omni_data, trange=None, keys=None, f=None, ax=None):
    """Plot a dynamic multi-panel OMNI overview from preloaded data."""
    if keys is None:
        keys = ["BZ_GSM", "Pressure", "Beta", "Mgs_mach_num", "Vx", "proton_density"]
    keys = list(keys)
    if not keys:
        raise ValueError("keys must contain at least one OMNI variable")
    specs = [OMNI_PANELS.get(key, PanelSpec((key,))) for key in keys]
    return plot_panels(omni_data, specs, trange=trange, f=f, ax=ax)


@tools.pass_axes()
def dst(dst_data, trange=None, f=None, ax=None):
    return _plot_panel(
        dst_data,
        PanelSpec(("dst",), zero=True),
        trange=trange,
        f=f,
        ax=ax,
    )

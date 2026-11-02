"""MMS particle panels built from declarative PanelSpecs."""

from plots import tools
from plots.panels import PanelSpec, plot_panel

__all__ = ["fpi", "hpca", "overview"]


_FPI_PANELS = (
    PanelSpec(("numberdensity_fast",), yscale="log"),
    PanelSpec(("tempperp_fast", "temppara_fast"), yscale="log", legend=True),
    PanelSpec(("bulkv_gse_fast",)),
)

_HPCA_PANELS = (
    PanelSpec(
        ("hplus_number_density", "heplus_number_density", "oplus_number_density"),
        yscale="log",
        legend=True,
    ),
    PanelSpec(
        ("hplus_tperp", "heplus_tperp", "oplus_tperp"),
        yscale="log",
        legend=True,
    ),
    PanelSpec(("hplus_ion_bulk_velocity_GSM",), legend=True),
    PanelSpec(("heplus_ion_bulk_velocity_GSM",), legend=True),
    PanelSpec(("oplus_ion_bulk_velocity_GSM",), legend=True),
)


def _plot_product(data, specs, f, ax):
    if len(ax) != len(specs):
        raise ValueError(f"expected {len(specs)} axes, got {len(ax)}")
    for axis, spec in zip(ax, specs):
        plot_panel(data, spec, f=f, ax=axis)
    return f, ax


@tools.pass_axes(subplot_spec=(3, 1, True, False))
def fpi(fpi, f=None, ax=None):
    return _plot_product(fpi, _FPI_PANELS, f, ax)


@tools.pass_axes(subplot_spec=(5, 1, True, False))
def hpca(hpca, f=None, ax=None):
    return _plot_product(hpca, _HPCA_PANELS, f, ax)


def overview(fpi_data=None, hpca_data=None, *, f=None, ax=None):
    """Compose whichever MMS particle products are available.

    FPI contributes three panels and HPCA five; passing both stacks all eight
    on a shared UTC axis. Supplied axes are reused, so callers can embed this
    in a larger report.
    """
    if fpi_data is None and hpca_data is None:
        raise ValueError("at least one of fpi_data or hpca_data is required")
    panel_count = (3 if fpi_data is not None else 0) + (5 if hpca_data is not None else 0)
    f, axes = tools.axes_for_columns(panel_count, f=f, ax=ax, sharex=True)
    index = 0
    if fpi_data is not None:
        fpi(fpi_data, f=f, ax=axes[index:index + 3])
        index += 3
    if hpca_data is not None:
        hpca(hpca_data, f=f, ax=axes[index:index + 5])
    axes[-1].set_xlabel("Time (UTC)")
    f.autofmt_xdate()
    return f, axes

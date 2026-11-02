"""
THEMIS product plots: multi-panel plots, energy spectra, and
magnetic-field / plasma time-series for THEMIS probes.

Uses shared plotting plumbing from ``plots.tools`` plus THEMIS-specific panel
helpers.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ── Decorator / plumbing from plots.tools ───────────────────────────
from libs import plasma
from plots import tools
from plots.panels import PanelSpec, plot_panel
from libs.spectrum import normalize_flux_by_energy_band
from .timeseries import dataset_variable, tplot as _tplot
from common.time import Time

__all__ = [
    "esa_spectrum",
    "esa_sst_spectra",
    "flux",
    "get_ma_beta",
    "overview",
    "spectrum",
    "sst_spectrum",
    "themis",
    "time_averaged_particle_moments",
    "time_averaged_spectra",
    "velocity_parallel_perp",
]


# ---------------------------------------------------------------------------
# Analysis helpers
# These operate on already-loaded frames. New report code should calculate
# statistics outside the plotting layer and pass the results to plot functions.
# ---------------------------------------------------------------------------


def get_ma_beta(fgm, esa, mom, prefix=None, thres_ma=50, thres_bt=50):
    """Return Mach number and beta on the FGM time grid.

    Invalid source samples are removed before interpolation. Values outside
    the source time ranges are NaN rather than extrapolated.
    ``prefix`` accepts either ``"thb"`` or the legacy ``"thb_"`` form.
    Pass ``None`` or ``np.inf`` for a threshold to disable display clipping.
    """
    _, fgm_B = dataset_variable(fgm, 'fgs_btotal', prefix)
    _, esa_n = dataset_variable(esa, 'peir_density', prefix)
    _, esa_T = dataset_variable(esa, 'peir_avgtemp', prefix)
    _, mom_v = dataset_variable(mom, 'peim_velocity_gsm', prefix)
    # Use the high-cadence magnetic-field time grid. Source NaNs must be
    # removed before interpolation; scipy.interp1d otherwise propagates them
    # across neighboring intervals and creates artificial gaps.
    t_eval = np.asarray(Time(fgm_B['x']).to_timestamp(), dtype=float)

    def interpolate_valid(variable, values):
        source_t = np.asarray(Time(variable['x']).to_timestamp(), dtype=float)
        source_y = np.asarray(values, dtype=float)
        valid = np.isfinite(source_t) & np.isfinite(source_y)
        if valid.sum() < 2:
            return np.full(t_eval.shape, np.nan, dtype=float)
        source_t, source_y = source_t[valid], source_y[valid]
        order = np.argsort(source_t)
        source_t, source_y = source_t[order], source_y[order]
        source_t, unique = np.unique(source_t, return_index=True)
        source_y = source_y[unique]
        return np.interp(t_eval, source_t, source_y, left=np.nan, right=np.nan)

    f_n = interpolate_valid(esa_n, esa_n['y'])
    f_T = interpolate_valid(esa_T, esa_T['y'])
    f_v = interpolate_valid(mom_v, np.linalg.norm(mom_v['y'], axis=1))
    ma = plasma.alfven_mach(fgm_B['y'], f_n, f_v)
    bt = plasma.plasma_beta(fgm_B['y'], f_n, f_T)

    ma[ma < 0] = np.nan
    if thres_ma is not None and np.isfinite(thres_ma):
        ma[ma > thres_ma] = np.nan
    bt[bt < 0] = np.nan
    if thres_bt is not None and np.isfinite(thres_bt):
        bt[bt > thres_bt] = np.nan
    return t_eval, ma, bt


def velocity_parallel_perp(fgm, mom, prefix=None, max_distance_seconds=60.0):
    """Project ion bulk velocity parallel/perpendicular to the local field.

    The FGM vector is interpolated onto the MOM velocity time grid. Parallel
    velocity is signed along ``B``; perpendicular velocity is the non-negative
    magnitude of the remaining vector component. Samples farther than
    ``max_distance_seconds`` from a valid magnetic-field record are left NaN.
    """
    _, magnetic = dataset_variable(fgm, 'fgs_gsm', prefix)
    _, velocity = dataset_variable(mom, 'peim_velocity_gsm', prefix)
    magnetic_time = np.asarray(Time(magnetic['x']).to_timestamp(), dtype=float)
    velocity_time = np.asarray(Time(velocity['x']).to_timestamp(), dtype=float)
    magnetic_value = np.asarray(magnetic['y'], dtype=float)
    velocity_value = np.asarray(velocity['y'], dtype=float)

    valid_b = np.isfinite(magnetic_time) & np.all(np.isfinite(magnetic_value), axis=1)
    if valid_b.sum() < 2:
        empty = np.full(len(velocity_time), np.nan, dtype=float)
        return pd.DatetimeIndex(velocity['x']), empty, empty.copy()
    source_time = magnetic_time[valid_b]
    source_b = magnetic_value[valid_b]
    order = np.argsort(source_time)
    source_time, source_b = source_time[order], source_b[order]
    source_time, unique = np.unique(source_time, return_index=True)
    source_b = source_b[unique]

    b_interp = np.column_stack([
        np.interp(velocity_time, source_time, source_b[:, component], left=np.nan, right=np.nan)
        for component in range(3)
    ])
    nearest = np.searchsorted(source_time, velocity_time)
    right = np.clip(nearest, 0, len(source_time) - 1)
    left = np.clip(nearest - 1, 0, len(source_time) - 1)
    distance = np.minimum(
        np.abs(velocity_time - source_time[left]),
        np.abs(velocity_time - source_time[right]),
    )
    b_interp[distance > max_distance_seconds] = np.nan

    b_norm = np.linalg.norm(b_interp, axis=1)
    valid = (
        np.isfinite(b_norm) & (b_norm > 0)
        & np.all(np.isfinite(velocity_value), axis=1)
    )
    b_hat = np.full_like(b_interp, np.nan, dtype=float)
    b_hat[valid] = b_interp[valid] / b_norm[valid, None]
    v_parallel = np.sum(velocity_value * b_hat, axis=1)
    v_squared = np.sum(velocity_value * velocity_value, axis=1)
    v_perp = np.sqrt(np.clip(v_squared - v_parallel * v_parallel, 0.0, None))
    v_parallel[~valid] = np.nan
    v_perp[~valid] = np.nan
    return pd.DatetimeIndex(velocity['x']), v_parallel, v_perp


def _spectrum_arrays(product, key, trange=None, prefix=None, *, en=True):
    """Validate and select a THEMIS spectrum.

    THEMIS L2 ``*_en_eflux`` variables contain differential energy flux.
    With ``en=False`` they are converted channel-by-channel to differential
    particle flux using ``j(E) = Eflux(E) / E``.
    """
    if not isinstance(en, (bool, np.bool_)):
        raise TypeError(f"en must be bool, got {type(en).__name__}")
    resolved, dat = dataset_variable(product, key, prefix=prefix)
    try:
        t = pd.DatetimeIndex(pd.to_datetime(dat["x"]))
        flux = np.asanyarray(dat["y"])
        energy = np.asanyarray(dat["v"])
    except KeyError as exc:
        raise ValueError(f"Spectrum {resolved!r} must contain x, y, and v arrays") from exc
    if flux.ndim != 2 or flux.shape[0] != len(t):
        raise ValueError(f"Spectrum {resolved!r}: y must have shape (time, energy), got {flux.shape}")
    if energy.ndim == 1:
        if energy.shape[0] != flux.shape[1]:
            raise ValueError(f"Spectrum {resolved!r}: v length {energy.shape[0]} does not match y energy dimension {flux.shape[1]}")
        energy = np.broadcast_to(energy, flux.shape)
    elif energy.ndim == 2:
        if energy.shape != flux.shape:
            raise ValueError(f"Spectrum {resolved!r}: v shape {energy.shape} does not match y shape {flux.shape}")
    else:
        raise ValueError(f"Spectrum {resolved!r}: v must be 1-D or 2-D, got {energy.shape}")
    mask = np.ones(len(t), dtype=bool)
    if trange is not None:
        if len(trange) != 2:
            raise ValueError("trange must contain exactly [start, end]")
        start, end = pd.to_datetime(trange[0]), pd.to_datetime(trange[1])
        mask = (t >= start) & (t <= end)
    if not mask.any():
        raise ValueError(f"Spectrum {resolved!r}: no samples in the requested time range")
    t, energy, flux = t[mask], energy[mask], flux[mask]
    flux = np.ma.filled(flux.astype(float), np.nan)
    energy = np.ma.filled(energy.astype(float), np.nan)
    if not en:
        flux = np.divide(
            flux,
            energy,
            out=np.full(flux.shape, np.nan, dtype=float),
            where=np.isfinite(energy) & (energy > 0.0),
        )
    if not np.isfinite(flux[flux > 0]).any():
        quantity = "energy-flux" if en else "particle-flux"
        raise ValueError(f"Spectrum {resolved!r}: no positive finite {quantity} values to plot")
    return resolved, t, energy, flux


@tools.pass_axes()
def spectrum(product, key, trange=None, prefix=None, en=True, f=None, ax=None, **kwargs):
    """Draw one energy spectrum and return its ``QuadMesh`` (no colorbar).

    ``en=True`` (default) plots differential energy flux. ``en=False`` plots
    differential particle flux. Positive raw flux values are displayed with a
    logarithmic normalization by default, so attached colorbars retain the
    original physical values rather than their base-10 exponents.
    """
    _, times, energy, flux = _spectrum_arrays(
        product, key, trange, prefix, en=en
    )
    # pcolormesh expects its first matrix axis to be y (energy) and its
    # second to be x (time), hence the transpose from the data's (time, bin)
    # storage convention.
    safe_energy = pd.DataFrame(energy).ffill().bfill().ffill(axis=1).bfill(axis=1).to_numpy().T
    z = np.where(np.isfinite(flux) & (flux > 0), flux, np.nan).T
    time_grid = np.broadcast_to(times.to_numpy()[None, :], safe_energy.shape)
    plot_kwargs = dict(kwargs)
    plot_kwargs.setdefault("norm", "log")
    return ax.pcolormesh(time_grid, safe_energy, z, shading="auto", rasterized=True, **plot_kwargs)


# ── Multi-panel THEMIS plots ───────────────────────────────────────────
_ESA_DENSITY = PanelSpec(("peir_density",), yscale="log")
_ESA_TEMP = PanelSpec(("peir_avgtemp",), yscale="log")
_FGM_BFIELD = PanelSpec(
    ("fgs_btotal", "fgs_gsm"),
    ylabel=r"$B_{GSM}$ (nT)",
    zero=True,
    legend=True,
    key_kwargs=({"color": "k"}, {}),
)
_MOM_VELOCITY = PanelSpec(
    ("peim_velocity_gsm",),
    ylabel=r"$V_{GSM}$ (km/s)",
    zero=True,
    legend=True,
)


def _style_spectrum_axis(ax, mesh, ylim, dir, *, en=True):
    # Keep the colorbar beside its own spectrum so it does not consume the
    # vertical space between the SST and ESA panels.
    cbar = ax.get_figure().colorbar(mesh, ax=ax, location=dir, pad=0.02)
    cbar.set_label(
        r"$Eflux$" if en else r"$j$"
    )
    cbar.ax.tick_params(labelsize=10)
    ax.set_yscale('log')
    ax.set_ylim(*ylim)
    ax.set_ylabel(r"$E \mathrm{(eV)}$",rotation=0,labelpad=30)

#recursive funcition, repeat args, how to eliminiate it
@tools.pass_axes()
def esa_spectrum(esa, ax=None, trange=None, *, key='peir_en_eflux', dir='top',
                 prefix=None, en=True, f=None, **kwargs):
    """Draw one ESA spectrum, including its single colorbar."""
    im = spectrum(
        esa, key, trange=trange, prefix=prefix, en=en, f=f, ax=ax, **kwargs
    )
    _style_spectrum_axis(ax, im, (10, 25000), dir=dir, en=en)
    return f, ax, im

@tools.pass_axes()
def sst_spectrum(sst, ax=None, trange=None, *, key='psif_en_eflux', dir='top',
                 prefix=None, en=True, f=None, **kwargs):
    """Draw one SST spectrum, including its single colorbar."""
    im = spectrum(
        sst, key, trange=trange, prefix=prefix, en=en, f=f, ax=ax, **kwargs
    )
    _style_spectrum_axis(ax, im, (30000, 200000), dir=dir, en=en)
    return f, ax, im

@tools.pass_axes((2, 1, True, False))
def esa_sst_spectra(esa, sst, trange=None, *, esa_key='peir_en_eflux',
                    sst_key='psif_en_eflux', dir='top',prefix=None,
                    shared_color_scale=False, en=True, f=None, ax=None,
                    **mesh_kwargs):
    """Plot SST (top) and ESA (bottom) spectra on separate panels.

    The instruments retain their own time and energy grids: this function
    never interpolates or concatenates them. Pass external axes as
    ``ax=(sst_ax, esa_ax)``; omitted axes are created by the decorator.
    """
    axes = np.asarray(ax, dtype=object).ravel().tolist()
    if len(axes) != 2:
        raise ValueError('ax must contain exactly (sst_ax, esa_ax)')
    sst_ax, esa_ax = axes
    if sst_ax.get_figure() is not f or esa_ax.get_figure() is not f:
        raise ValueError('ax and f must belong to the same figure')

    # Validate both inputs before adding either mesh, so malformed data does
    # not leave a half-rendered figure behind.
    _, _, _, esa_flux = _spectrum_arrays(
        esa, esa_key, trange, prefix, en=en
    )
    _, _, _, sst_flux = _spectrum_arrays(
        sst, sst_key, trange, prefix, en=en
    )
    plot_kwargs = dict(mesh_kwargs)
    if shared_color_scale:
        positive = np.concatenate((esa_flux[esa_flux > 0], sst_flux[sst_flux > 0]))
        limits = (np.nanmin(positive), np.nanmax(positive))
        plot_kwargs.setdefault('vmin', limits[0])
        plot_kwargs.setdefault('vmax', limits[1])

    sst_mesh = spectrum(
        sst, sst_key, trange=trange, prefix=prefix, en=en,
        f=f, ax=sst_ax, **plot_kwargs
    )
    esa_mesh = spectrum(
        esa, esa_key, trange=trange, prefix=prefix, en=en,
        f=f, ax=esa_ax, **plot_kwargs
    )
    _style_spectrum_axis(
        sst_ax, sst_mesh, (30000, 200000), dir=dir, en=en
    )
    _style_spectrum_axis(
        esa_ax, esa_mesh, (10, 25000), dir=dir, en=en
    )
    esa_ax.set_xlabel('Time (UTC)')
    return f, (sst_ax, esa_ax), (sst_mesh, esa_mesh)

#fgm
#B field gsm and btotal spin resolution
def themis(esa, mom, fgm, panels=5, trange=None, en=True, f=None, ax=None):
    if panels not in (3, 5, 7):
        raise ValueError(f"themis() supports panels in (3, 5, 7), got {panels!r}")
    if ax is None:
        f, axes = plt.subplots(panels, 1, sharex=True)
        ax = axes.flatten().tolist()
    else:
        ax = np.asarray(ax, dtype=object).ravel().tolist()
        if len(ax) != panels:
            raise ValueError(f'themis() requires exactly {panels} axes, got {len(ax)}')
        if f is None:
            f = ax[0].get_figure()
    if panels == 3:
        plot_panel(esa, _ESA_DENSITY, trange=trange, f=f, ax=ax[0])
        axb = ax[0].twinx()
        _tplot(esa, 'peir_avgtemp', trange=trange, f=f, ax=axb, color="r")
        ax[0].set_yscale("log")
        axb.set_yscale("log")
        ax[0].axhline(0.03, ls="--")
        plot_panel(mom, _MOM_VELOCITY, trange=trange, f=f, ax=ax[1])
        plot_panel(fgm, _FGM_BFIELD, trange=trange, f=f, ax=ax[2])
    elif panels == 5:
        plot_panel(esa, _ESA_DENSITY, trange=trange, f=f, ax=ax[0])
        plot_panel(esa, _ESA_TEMP, trange=trange, f=f, ax=ax[1])
        for i in range(2):
            ax[i].set_yscale("log")
        ax[0].axhline(0.03, ls="--")
        esa_spectrum(esa, ax[2], trange=trange, en=en, f=f)
        plot_panel(mom, _MOM_VELOCITY, trange=trange, f=f, ax=ax[3])
        plot_panel(fgm, _FGM_BFIELD, trange=trange, f=f, ax=ax[4])
    elif panels == 7:
        plot_panel(esa, _ESA_DENSITY, trange=trange, f=f, ax=ax[0])
        plot_panel(esa, _ESA_TEMP, trange=trange, f=f, ax=ax[1])
        for i in range(2):
            ax[i].set_yscale("log")
        ax[0].axhline(0.03, ls="--")
        esa_spectrum(esa, ax[2], trange=trange, en=en, f=f)
        plot_panel(mom, _MOM_VELOCITY, trange=trange, f=f, ax=ax[3])
        plot_panel(fgm, _FGM_BFIELD, trange=trange, f=f, ax=ax[4])
        t, ma, beta = get_ma_beta(fgm, esa, mom)
        ax[5].plot(pd.to_datetime(t, unit="s"), ma)
        ax[6].plot(pd.to_datetime(t, unit="s"), beta)
        for i in (5, 6):
            ax[i].set_yscale("log")
        ax[6].set_xlabel("Time (UTC)")
    return f, ax


def overview(esa=None, mom=None, fgm=None, *, trange=None, en=True, f=None, ax=None):
    """Compose the THEMIS products that are actually available.

    ``themis()`` is the fixed 3/5/7-panel layout and requires all three
    products. Quick-look callers often have only one or two products; this
    helper builds a vertical stack from whichever of ``esa``, ``mom`` and
    ``fgm`` are supplied, preserving each product's native time axis.
    """
    if esa is None and mom is None and fgm is None:
        raise ValueError("at least one of esa, mom or fgm is required")

    panel_count = (3 if esa is not None else 0)
    panel_count += 1 if mom is not None else 0
    panel_count += 1 if fgm is not None else 0
    f, axes = tools.axes_for_columns(panel_count, f=f, ax=ax, sharex=True)
    index = 0
    if esa is not None:
        plot_panel(esa, _ESA_DENSITY, trange=trange, f=f, ax=axes[index])
        axes[index].set_yscale("log")
        index += 1
        plot_panel(esa, _ESA_TEMP, trange=trange, f=f, ax=axes[index])
        axes[index].set_yscale("log")
        index += 1
        esa_spectrum(esa, axes[index], trange=trange, en=en, f=f)
        index += 1
    if mom is not None:
        plot_panel(mom, _MOM_VELOCITY, trange=trange, f=f, ax=axes[index])
        index += 1
    if fgm is not None:
        plot_panel(fgm, _FGM_BFIELD, trange=trange, f=f, ax=axes[index])
        index += 1

    axes[-1].set_xlabel("Time (UTC)")
    f.autofmt_xdate()
    return f, axes


# ---------------------------------------------------------------------------
# Spectrum and moment analysis retained below for existing THEMIS workflows.
# ---------------------------------------------------------------------------
def time_averaged_spectra(esa, sst, tspan=None, *, esa_key='peir_en_eflux',
                          sst_key='psif_en_eflux', prefix=None, en=True,
                          normalized=False, E_norm=None):
    """Return one cleaned, energy-sorted ESA+SST ``(energy, flux)`` pair.

    ``en=True`` returns differential energy flux. ``en=False`` first converts
    every time/energy sample to differential particle flux and then averages,
    so moving energy-bin centers are handled correctly.

    Energy channels for which either time-averaged value is non-finite are
    removed as a pair before the ESA and SST points are concatenated and
    sorted. No interpolation is performed across the instrument gap.

    Set ``normalized=True`` to divide the combined spectrum by its mean flux
    in ``E_norm``. THEMIS Energy values are in eV, so use
    ``E_norm=(10000, 15000)`` for 10--15 keV.
    """
    if not isinstance(normalized, (bool, np.bool_)):
        raise TypeError(f"normalized must be bool, got {type(normalized).__name__}")
    if normalized and E_norm is None:
        raise ValueError("E_norm is required when normalized=True.")
    if not normalized and E_norm is not None:
        raise ValueError("E_norm is only valid when normalized=True.")

    def average(product, key):
        _, times, energy, flux = _spectrum_arrays(
            product, key, None, prefix, en=en
        )
        if tspan is not None:
            start, end = pd.to_datetime(tspan[0]), pd.to_datetime(tspan[1])
            selected = (times >= start) & (times <= end)
            if not selected.any():
                return None
            energy, flux = energy[selected], flux[selected]
        energy = np.ma.filled(np.asanyarray(energy, dtype=float), np.nan)
        flux = np.ma.filled(np.asanyarray(flux, dtype=float), np.nan)
        energy_count = np.isfinite(energy).sum(axis=0)
        flux_count = np.isfinite(flux).sum(axis=0)
        energy_mean = np.divide(
            np.nansum(energy, axis=0),
            energy_count,
            out=np.full(energy.shape[1], np.nan),
            where=energy_count > 0,
        )
        flux_mean = np.divide(
            np.nansum(flux, axis=0),
            flux_count,
            out=np.full(flux.shape[1], np.nan),
            where=flux_count > 0,
        )
        valid = np.isfinite(energy_mean) & np.isfinite(flux_mean)
        return energy_mean[valid], flux_mean[valid]

    averages = [average(esa, esa_key), average(sst, sst_key)]
    available = [
        spectrum for spectrum in averages
        if spectrum is not None and spectrum[0].size
    ]
    if available:
        energy = np.concatenate([spectrum[0] for spectrum in available])
        flux = np.concatenate([spectrum[1] for spectrum in available])
        order = np.argsort(energy)
        energy, flux = energy[order], flux[order]
        if normalized:
            flux, _ = normalize_flux_by_energy_band(
                energy, flux, E_norm=E_norm
            )
        return energy, flux
    return None


def _time_averaged_product_moments(
    product,
    tspan=None,
    *,
    density_key='peir_density',
    velocity_key='peir_velocity_gsm',
    temperature_key='peir_avgtemp',
    prefix=None,
):
    """Return moments for one time-averaged particle product.

    Density is time averaged and velocity is density weighted as
    ``V = <n V> / <n>``. The scalar temperature also includes the spread
    caused by changes of bulk velocity inside the selected interval.
    """

    def samples(key, components=None):
        resolved, variable = dataset_variable(product, key, prefix=prefix)
        try:
            times = pd.DatetimeIndex(pd.to_datetime(variable['x']))
            values = np.ma.filled(
                np.asanyarray(variable['y'], dtype=float), np.nan
            )
        except KeyError as exc:
            raise ValueError(
                f"Moment {resolved!r} must contain x and y arrays"
            ) from exc
        expected = (len(times),) if components is None else (
            len(times), components
        )
        if values.shape != expected:
            raise ValueError(
                f"Moment {resolved!r}: y must have shape {expected}, "
                f"got {values.shape}"
            )
        order = np.argsort(times.view('i8'))
        return times[order], values[order]

    density_time, density = samples(density_key)
    velocity_time, velocity = samples(velocity_key, components=3)
    temperature_time, temperature = samples(temperature_key)

    clip_duration_to_span = tspan is not None
    if tspan is None:
        start = max(
            density_time.min(), velocity_time.min(), temperature_time.min()
        )
        stop = min(
            density_time.max(), velocity_time.max(), temperature_time.max()
        )
    else:
        if len(tspan) != 2:
            raise ValueError("tspan must contain exactly [start, end]")
        start, stop = pd.to_datetime(tspan[0]), pd.to_datetime(tspan[1])
    selected = (velocity_time >= start) & (velocity_time <= stop)
    if not selected.any():
        raise ValueError("No velocity samples in the requested time range")

    evaluation_time = velocity_time[selected]
    velocity = velocity[selected]
    origin_ns = min(
        density_time[0].value,
        velocity_time[0].value,
        temperature_time[0].value,
    )

    def interpolate(source_time, source_value):
        source_seconds = (
            source_time.view('i8').astype(float) - origin_ns
        ) * 1e-9
        evaluation_seconds = (
            evaluation_time.view('i8').astype(float) - origin_ns
        ) * 1e-9
        valid = np.isfinite(source_seconds) & np.isfinite(source_value)
        if valid.sum() < 2:
            return np.full(len(evaluation_time), np.nan)
        source_seconds = source_seconds[valid]
        source_value = source_value[valid]
        source_seconds, unique = np.unique(
            source_seconds, return_index=True
        )
        source_value = source_value[unique]
        return np.interp(
            evaluation_seconds,
            source_seconds,
            source_value,
            left=np.nan,
            right=np.nan,
        )

    density = interpolate(density_time, density)
    temperature = interpolate(temperature_time, temperature)
    valid = (
        np.isfinite(density)
        & (density > 0.0)
        & np.isfinite(temperature)
        & (temperature > 0.0)
        & np.all(np.isfinite(velocity), axis=1)
    )
    if valid.sum() == 0:
        raise ValueError("No positive finite particle moments in time range")

    evaluation_time = evaluation_time[valid]
    density = density[valid]
    temperature = temperature[valid]
    velocity = velocity[valid]

    if len(evaluation_time) == 1:
        duration = np.ones(1)
    else:
        seconds = evaluation_time.view('i8').astype(float) * 1e-9
        edges = np.empty(len(seconds) + 1)
        edges[1:-1] = 0.5 * (seconds[:-1] + seconds[1:])
        edges[0] = seconds[0] - 0.5 * (seconds[1] - seconds[0])
        edges[-1] = seconds[-1] + 0.5 * (
            seconds[-1] - seconds[-2]
        )
        if clip_duration_to_span:
            edges[0] = max(
                pd.Timestamp(start).value * 1e-9, edges[0]
            )
            edges[-1] = min(
                pd.Timestamp(stop).value * 1e-9, edges[-1]
            )
        duration = np.maximum(np.diff(edges), 0.0)
        if not np.any(duration > 0.0):
            duration = np.ones(len(seconds))

    density_mean = np.average(density, weights=duration)
    density_weight = duration * density
    velocity_mean = np.average(
        velocity, axis=0, weights=density_weight
    )
    intrinsic_temperature = np.average(
        temperature, weights=density_weight
    )

    from common.constants import mass_H_kg, e_charge
    velocity_variance_m2s2 = np.average(
        np.sum((velocity - velocity_mean)**2, axis=1) * 1e6,
        weights=density_weight,
    )
    mixing_temperature = (
        mass_H_kg * velocity_variance_m2s2 / (3.0 * e_charge)
    )

    return {
        'density_cm3': float(density_mean),
        'velocity_gsm_kms': np.asarray(velocity_mean, dtype=float),
        'speed_kms': float(np.linalg.norm(velocity_mean)),
        'temperature_eV': float(
            intrinsic_temperature + mixing_temperature
        ),
        'intrinsic_temperature_eV': float(intrinsic_temperature),
        'temporal_mixing_temperature_eV': float(mixing_temperature),
        'sample_count': int(valid.sum()),
    }


def time_averaged_particle_moments(
    esa,
    sst=None,
    tspan=None,
    *,
    esa_density_key='peir_density',
    esa_velocity_key='peir_velocity_gsm',
    esa_temperature_key='peir_avgtemp',
    sst_density_key='psif_density',
    sst_velocity_key='psif_velocity_gsm',
    sst_temperature_key='psif_avgtemp',
    prefix=None,
):
    """Return moments matching an ESA-only or combined ESA+SST spectrum.

    ESA and SST cover different ion-energy ranges. When ``sst`` is supplied,
    their partial density, momentum, and scalar pressure moments are added.
    The between-instrument velocity separation is included in the combined
    scalar temperature.
    """
    esa_moments = _time_averaged_product_moments(
        esa,
        tspan,
        density_key=esa_density_key,
        velocity_key=esa_velocity_key,
        temperature_key=esa_temperature_key,
        prefix=prefix,
    )
    if sst is None:
        result = dict(esa_moments)
        result['instruments'] = {'esa': esa_moments}
        result['instrument_mixing_temperature_eV'] = 0.0
        return result

    sst_moments = _time_averaged_product_moments(
        sst,
        tspan,
        density_key=sst_density_key,
        velocity_key=sst_velocity_key,
        temperature_key=sst_temperature_key,
        prefix=prefix,
    )
    parts = (esa_moments, sst_moments)
    densities = np.array([part['density_cm3'] for part in parts])
    velocities = np.vstack([
        part['velocity_gsm_kms'] for part in parts
    ])
    temperatures = np.array([
        part['temperature_eV'] for part in parts
    ])
    total_density = densities.sum()
    if not np.isfinite(total_density) or total_density <= 0.0:
        raise ValueError("Combined ESA+SST density is not positive")
    total_velocity = np.sum(
        densities[:, None] * velocities, axis=0
    ) / total_density

    from common.constants import mass_H_kg, e_charge
    between_variance_m2s2 = np.sum(
        densities
        * np.sum((velocities - total_velocity)**2, axis=1)
    ) / total_density * 1e6
    instrument_mixing_temperature = (
        mass_H_kg * between_variance_m2s2 / (3.0 * e_charge)
    )
    within_temperature = np.sum(
        densities * temperatures
    ) / total_density
    intrinsic_temperature = np.sum(
        densities
        * np.array([
            part['intrinsic_temperature_eV'] for part in parts
        ])
    ) / total_density
    temporal_mixing_temperature = np.sum(
        densities
        * np.array([
            part['temporal_mixing_temperature_eV'] for part in parts
        ])
    ) / total_density

    return {
        'density_cm3': float(total_density),
        'velocity_gsm_kms': np.asarray(total_velocity, dtype=float),
        'speed_kms': float(np.linalg.norm(total_velocity)),
        'temperature_eV': float(
            within_temperature + instrument_mixing_temperature
        ),
        'intrinsic_temperature_eV': float(intrinsic_temperature),
        'temporal_mixing_temperature_eV': float(
            temporal_mixing_temperature
        ),
        'instrument_mixing_temperature_eV': float(
            instrument_mixing_temperature
        ),
        'sample_count': int(sum(
            part['sample_count'] for part in parts
        )),
        'instruments': {
            'esa': esa_moments,
            'sst': sst_moments,
        },
    }


def flux(esa, sst, prob=None, tspan=None, plot=False, f=None, ax=None, *,
          esa_key='peir_en_eflux', sst_key='psif_en_eflux', prefix=None,
          en=True, normalized=False, E_norm=None):
    """Compute one cleaned ESA+SST average spectrum and optionally plot it.

    The historical name is retained for compatibility. Set ``en=False`` for
    differential particle flux. Set ``normalized=True`` with ``E_norm`` in
    eV to compare spectral shapes instead of absolute amplitudes.
    """
    if prob is not None:
        prefix = f'th{prob}'
    average = time_averaged_spectra(esa, sst, tspan, esa_key=esa_key,
                                    sst_key=sst_key, prefix=prefix, en=en,
                                    normalized=normalized, E_norm=E_norm)
    if plot:
        if ax is None:
            f, ax = plt.subplots() if f is None else (f, f.subplots())
        elif f is None:
            f = ax.get_figure()
        if average is not None:
            energy, flux = average
            ax.plot(
                energy, flux, 'k+', linestyle='none', markersize=5,
                label=(
                    (
                        'Normalized ESA+SST energy flux'
                        if en else 'Normalized ESA+SST particle flux'
                    ) if normalized
                    else (
                        'Time-averaged ESA+SST energy flux'
                        if en else 'Time-averaged ESA+SST particle flux'
                    )
                ),
            )
        if normalized:
            ylabel = r'Normalized energy flux' if en else r'Normalized particle flux'
        elif en:
            ylabel = r'Differential energy flux $(eV\,cm^{-2}\,s^{-1}\,sr^{-1}\,eV^{-1})$'
        else:
            ylabel = r'Differential particle flux $(cm^{-2}\,s^{-1}\,sr^{-1}\,eV^{-1})$'
        ax.set(
            xscale='log', yscale='log', xlabel='Energy (eV)',
            ylabel=ylabel,
        )
        ax.grid(True, which='both', linestyle='--', alpha=0.7)
        ax.legend()
        f.tight_layout()
        return f, ax
    else:
        return average

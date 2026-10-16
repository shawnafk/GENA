"""Post-processing statistics for GENA particle rows.

Functions in this module operate on ordinary pandas DataFrames and return
derived statistics, spectra, maps, or fitted quantities. ParticleDataFrame keeps
methods that return DataFrame subsets or add standard event-level columns.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Callable

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from scipy import stats

from common.constants import e_charge, mass_H_kg
from libs.spectrum import normalize_flux_by_energy_band
from common.time import Time
from libs.plasma import sigma_HH


SCIME_ENERGY_BINS_KEV = np.geomspace(4.0, 25.0, 11)
SCIME_MIN_ENERGY_KEV = 8.0
SCIME_MAX_ENERGY_KEV = 25.0


@dataclass
class CountsNDResult:
    """N-dimensional raw and optionally weighted event counts."""

    dims: tuple[str, ...]
    centers: dict[str, np.ndarray]
    edges: dict[str, np.ndarray]
    collapsed_ranges: dict[str, tuple]
    raw_counts: np.ndarray
    raw_count_errors: np.ndarray
    weighted_counts: np.ndarray | None
    weighted_count_errors: np.ndarray | None
    input_size: int
    selected_size: int


def counts_nd(
    df: pd.DataFrame,
    axes: Mapping[str, Sequence],
    *,
    weights: ArrayLike | None = None,
    keepdims: Sequence[str] = (),
) -> CountsNDResult:
    """Count events after range selection on any number of DataFrame columns.

    An axis with exactly two edges is selected over ``[low, high)`` and then
    collapsed unless its name appears in ``keepdims``. An axis with more than
    two edges is retained in the result. Retained result dimensions follow the
    insertion order of ``axes``.
    ``weights`` must be an external one-dimensional array aligned positionally
    with ``df``. No exposure, bin-width, geometric-factor, or other flux
    normalization is applied.
    """
    weight_values = None if weights is None else np.asarray(weights, dtype=float)
    if weight_values is not None and (weight_values.ndim != 1 or len(weight_values) != len(df)):
        raise ValueError("weights must be one-dimensional and aligned with df.")
    retained_names = set(keepdims)
    unknown_retained = retained_names.difference(axes)
    if unknown_retained:
        names = ", ".join(sorted(unknown_retained))
        raise ValueError(f"keepdims contains axes not present in axes: {names}")
    selected = np.ones(len(df), dtype=bool)
    values_by_axis = {}
    bins_by_axis = {}
    dims = []
    centers = {}
    output_edges = {}
    collapsed_ranges: dict[str, tuple] = {}

    for name, explicit_bins in axes.items():
        if name == "Datetime":
            value_ns = Time(df[name]).datetime.asi8
            edge_ns = Time(explicit_bins).datetime.asi8
            origin_ns = edge_ns[0]
            values = (value_ns - origin_ns).astype(float) / 1.0e9
            bins = (edge_ns - origin_ns).astype(float) / 1.0e9
            shown_edges = edge_ns.astype("datetime64[ns]")
            axis_centers = (edge_ns[:-1] + (edge_ns[1:] - edge_ns[:-1]) // 2).astype("datetime64[ns]")
        else:
            values = np.asarray(df[name], dtype=float)
            bins = np.asarray(explicit_bins, dtype=float)
            shown_edges = bins
            axis_centers,_,_ = get_bin_grid(bins)

        selected &= (values >= bins[0]) & (values < bins[-1])
        if len(bins) == 2 and name not in retained_names:
            collapsed_ranges[name] = (shown_edges[0], shown_edges[-1])
        else:
            dims.append(name)
            values_by_axis[name] = values
            bins_by_axis[name] = bins
            centers[name] = axis_centers
            output_edges[name] = shown_edges

    selected_size = int(np.sum(selected))

    if dims:
        samples = np.column_stack([values_by_axis[name][selected] for name in dims])
        histogram_bins = [bins_by_axis[name] for name in dims]
        raw_counts_float, _ = np.histogramdd(samples, bins=histogram_bins)
        raw_counts = raw_counts_float.astype(np.int64)
        if weight_values is None:
            weighted_counts = None
            weighted_count_errors = None
        else:
            selected_weights = weight_values[selected]
            weighted_counts, _ = np.histogramdd(samples, bins=histogram_bins, weights=selected_weights)
            weighted_variance, _ = np.histogramdd(samples, bins=histogram_bins, weights=selected_weights**2)
            weighted_count_errors = np.sqrt(weighted_variance)
    else:
        raw_counts = np.asarray(selected_size, dtype=np.int64)
        if weight_values is None:
            weighted_counts = None
            weighted_count_errors = None
        else:
            selected_weights = weight_values[selected]
            weighted_counts = np.asarray(np.sum(selected_weights), dtype=float)
            weighted_count_errors = np.asarray(np.sqrt(np.sum(np.square(selected_weights))), dtype=float)

    return CountsNDResult(
        dims=tuple(dims),
        centers=centers,
        edges=output_edges,
        collapsed_ranges=collapsed_ranges,
        raw_counts=raw_counts,
        raw_count_errors=np.sqrt(raw_counts),
        weighted_counts=weighted_counts,
        weighted_count_errors=weighted_count_errors,
        input_size=len(df),
        selected_size=selected_size,
    )


def df_resample(df: pd.DataFrame, window: str = "5min"):
    return pd.DataFrame(df).resample(window, on="Datetime")

def ave(y, **kwargs):
    """Return the mean and its standard error.

    ``weights`` are treated as fixed relative/reliability weights.  The
    weighted SEM uses the Kish effective sample size and reduces to
    ``scipy.stats.sem`` when all weights are equal.
    """
    y = np.asarray(y)
    weights = kwargs.get("weights")

    if weights is None:
        # 均值和标准误必须使用同一组有限样本。
        y = y[np.isfinite(y)]
        if len(y) == 0:
            return np.nan, np.nan
        return np.mean(y), stats.sem(y)

    # 加权样本中，数值和对应权重任一无效时都排除该事件。
    valid = np.isfinite(y) & np.isfinite(weights)
    y = y[valid]
    w = weights[valid]
    if len(y) == 0:
        return np.nan, np.nan

    sum_w = np.sum(w)
    sum_w2 = np.sum(w**2)
    if sum_w <= 0 or sum_w2 == 0:
        return np.nan, np.nan

    mean_val = np.average(y, weights=w)

    # 加权总体方差；这里的权重只出现一次。
    population_variance = np.sum(w * (y - mean_val) ** 2) / sum_w

    # 1 - sum(w²)/sum(w)² 是有限样本的无偏修正项。
    correction = 1 - sum_w2 / sum_w**2
    if correction <= 0:
        return mean_val, np.nan

    # Kish 有效样本量：权重越不均匀，独立信息量越小。
    effective_size = sum_w**2 / sum_w2
    sample_variance = population_variance / correction
    sem_val = np.sqrt(sample_variance / effective_size)
    return mean_val, sem_val

def hist(y, bins, **kwargs):
    """Weighted histogram with Poisson counting uncertainty.

    Use ``density=False`` for counts that will be converted to physical flux.
    ``density=True`` returns a normalized probability density and already
    divides by both the total count and each bin width.
    """
    density = kwargs.pop("density", False)
    weights = kwargs.get("weights")

    # 始终先计算未归一化的计数
    w_bin, bin_edges = np.histogram(
        y,
        bins=bins,
        density=False,
        **kwargs,
    )

    if weights is None:
        # 未加权泊松方差：Var(N) = N
        err = np.sqrt(w_bin)
    else:
        # 加权泊松方差：Var(S) = sum(w_i²)
        w2_bin, _ = np.histogram( y, bins=bin_edges, weights=weights**2)
        err = np.sqrt(w2_bin)

    if density:
        widths = np.diff(bin_edges)

        # w_bin.sum() 只包含实际落入 bins 的事件或权重
        denominator = np.clip( w_bin.sum() * widths, 1e-30, None)

        w_bin = w_bin / denominator
        err = err / denominator

    centers, _, _ = get_bin_grid(bin_edges)
    return centers, w_bin, err


def hist2d(x, y, bins, **kwargs):
    """Weighted two-dimensional histogram with Poisson counting uncertainty.

    Use ``density=False`` before physical flux normalization. ``density=True``
    returns a probability density divided by the total count and bin area.
    """
    density = kwargs.pop("density", False)
    weights = kwargs.get("weights")

    counts, x_edges, y_edges = np.histogram2d(
        x,
        y,
        bins=bins,
        density=False,
        **kwargs,
    )

    if weights is None:
        err = np.sqrt(counts)
    else:
        counts2, _, _ = np.histogram2d(
            x,
            y,
            bins=[x_edges, y_edges],
            weights=weights**2,
        )
        err = np.sqrt(counts2)

    if density:
        areas = np.diff(x_edges)[:, None] * np.diff(y_edges)[None, :]
        denominator = np.clip(counts.sum() * areas, 1e-30, None)
        counts = counts / denominator
        err = err / denominator

    x_centers, _, _ = get_bin_grid(x_edges)
    y_centers, _, _ = get_bin_grid(y_edges)
    return x_centers, y_centers, counts, err

def window_ave(df: pd.DataFrame, key: str, window: str = "5min", **kwargs):
    """Return time-windowed averages and standard errors for one column."""
    weights = kwargs.get("weights")

    time_edges = []
    values = []
    errors = []
    for time, group in df_resample(df, window):
        if group.empty:
            continue
        group_kwargs = kwargs.copy()
        if weights is not None:
            group_kwargs["weights"] = weights[group.index]
        value, err = ave(group[key], **group_kwargs)
        time_edges.append(time)
        values.append(value)
        errors.append(err)
    return np.array(time_edges), np.array(values), np.array(errors)

def window_hist0d(df: pd.DataFrame, window: str = "5min", **kwargs):
    """Return event counts and Poisson errors in each time window."""
    weights = kwargs.get("weights")

    time_edges = []
    ret_hist = []
    errors = []
    for time, group in df_resample(df, window):
        if weights is None:
            count = len(group)
            err = np.sqrt(count)
        else:
            group_weights = weights[group.index]
            count = np.sum(group_weights)
            err = np.sqrt(np.sum(group_weights**2))
        ret_hist.append(count)
        errors.append(err)
        time_edges.append(time)
    return np.array(time_edges), np.array(ret_hist), np.array(errors)


def window_hist1d(df: pd.DataFrame, value_name: str, value_bins, window: str = "5min", **kwargs):
    """Return bin centers, histograms, and errors in each time window."""
    weights = kwargs.get("weights")

    time_edges = []
    centers = None
    histed = []
    errors = []
    for time, group in df_resample(df, window):
        group_kwargs = kwargs.copy()
        if weights is not None:
            group_kwargs["weights"] = weights[group.index]
        centers, counts, err = hist(group[value_name], value_bins, **group_kwargs)
        histed.append(counts)
        errors.append(err)
        time_edges.append(time)
    return np.array(time_edges), centers, np.array(histed), np.array(errors)


def window_hist2d(df: pd.DataFrame, value_names, value_bins, window: str = "5min", **kwargs):
    """Return axis centers, two-dimensional histograms, and errors by time."""
    weights = kwargs.get("weights")

    time_edges = []
    x_centers = None
    y_centers = None
    histed = []
    errors = []
    for time, group in df_resample(df, window):
        group_kwargs = kwargs.copy()
        if weights is not None:
            group_kwargs["weights"] = weights[group.index]
        x_centers, y_centers, counts, err = hist2d(
            group[value_names[0]],
            group[value_names[1]],
            value_bins,
            **group_kwargs,
        )
        histed.append(counts)
        errors.append(err)
        time_edges.append(time)
    return np.array(time_edges), x_centers, y_centers, np.array(histed), np.array(errors)


def get_bin_grid(bins, atol=1e-2):
    """Return centers, widths, and spacing type for explicit bin edges."""
    edges = np.asarray(bins, dtype=float)
    widths = np.diff(edges)
    if len(widths) == 0:
        raise ValueError("At least two bin edges are required.")

    all_linear = np.allclose(widths, widths[0], atol=atol)
    ratios = edges[1:] / edges[:-1] if np.all(edges > 0) else np.array([])
    all_log = len(ratios) > 1 and np.allclose(ratios, ratios[0], atol=atol)
    if all_linear:
        centers = (edges[:-1] + edges[1:]) / 2
        bin_type = "linear"
    elif all_log:
        centers = np.sqrt(edges[:-1] * edges[1:])
        bin_type = "log"
    else:
        centers = (edges[:-1] + edges[1:]) / 2
        bin_type = "irregular"
    return centers, widths, bin_type


def t_exposure(df):
    end = df.index[-1]
    start = df.index[0]
    return (df.loc[end,'Datetime']-df.loc[start,'Datetime']).total_seconds()



# -------------------------------------------------------------
# Counts-based flux and temperature inference
# -------------------------------------------------------------
def counts_to_flux(
    result: CountsNDResult,
    *,
    exposure: ArrayLike = 1.0,
    geometric_factor: float = 1.0,
    en: bool = False,
    normalized: bool = False,
    E_norm: Sequence | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert N-dimensional counts to differential particle flux.

    ``result`` must retain or collapse an ``Energy`` axis. Weighted counts are
    used when the caller supplied ``weights`` to :func:`counts_nd`; otherwise
    raw counts are used. Energy-bin width, exposure, and geometric factor are
    removed from counts and their matching statistical errors.

    Return ``(energy, flux, flux_error)``. Flux and its statistical error are
    divided by the same normalization. ``en=True`` multiplies both by the
    corresponding Energy centers. With ``normalized=True``, each retained
    Energy spectrum is divided by its mean returned flux over ``E_norm``
    (keV). The returned values are then dimensionless relative spectra.
    """
    if result.weighted_counts is None:
        counts = result.raw_counts
        count_errors = result.raw_count_errors
        if count_errors is None or np.shape(count_errors) != np.shape(counts):
            raise ValueError("raw counts and errors must both exist with matching shapes.")
    else:
        counts = result.weighted_counts
        count_errors = result.weighted_count_errors
        if count_errors is None or np.shape(count_errors) != np.shape(counts):
            raise ValueError("weighted counts and errors must both exist with matching shapes.")

    counts = np.asarray(counts, dtype=float)
    count_errors = np.asarray(count_errors, dtype=float)
    exposure_values = np.asarray(exposure, dtype=float)
    if np.any(~np.isfinite(exposure_values)) or np.any(exposure_values <= 0):
        raise ValueError("exposure must contain only positive finite values.")
    geometric_factor_value = np.asarray(geometric_factor, dtype=float)
    if geometric_factor_value.ndim != 0 or not np.isfinite(geometric_factor_value) or geometric_factor_value <= 0:
        raise ValueError("geometric_factor must be a positive finite scalar.")
    geometric_factor = float(geometric_factor_value)

    if "Energy" in result.dims:
        energy = np.asarray(result.centers["Energy"], dtype=float)
        energy_widths = np.diff(np.asarray(result.edges["Energy"], dtype=float))
        energy_shape = [1] * counts.ndim
        energy_shape[result.dims.index("Energy")] = len(energy)
        energy_widths = energy_widths.reshape(energy_shape)
        energy_factors = energy.reshape(energy_shape)
    elif "Energy" in result.collapsed_ranges:
        energy_edges = np.asarray(result.collapsed_ranges["Energy"], dtype=float)
        energy, energy_widths, _ = get_bin_grid(energy_edges)
        energy_widths = energy_widths[0]
        energy_factors = energy[0]
    else:
        raise ValueError("counts_to_flux requires an Energy axis or collapsed range.")

    try:
        exposure_values = np.broadcast_to(exposure_values, counts.shape)
    except ValueError as exc:
        raise ValueError("exposure is not broadcastable to the count shape.") from exc
    denominator = energy_widths * exposure_values * geometric_factor
    flux = counts / denominator
    flux_errors = count_errors / denominator
    if en:
        flux = flux * energy_factors
        flux_errors = flux_errors * energy_factors
    if not isinstance(normalized, (bool, np.bool_)):
        raise TypeError(f"normalized must be bool, got {type(normalized).__name__}")
    if normalized:
        if E_norm is None:
            raise ValueError("E_norm is required when normalized=True.")
        if "Energy" not in result.dims:
            raise ValueError("normalized=True requires a retained Energy axis.")
        flux, flux_errors = normalize_flux_by_energy_band(
            energy,
            flux,
            flux_errors,
            energy_axis=result.dims.index("Energy"),
            energy_edges=result.edges["Energy"],
            E_norm=E_norm,
        )
    elif E_norm is not None:
        raise ValueError("E_norm is only valid when normalized=True.")
    return energy, flux, flux_errors


def flux_e(
    df: pd.DataFrame,
    energy_bins: Sequence,
    *,
    exposure: ArrayLike = 1.0,
    geometric_factor: float = 1.0,
    weights: ArrayLike | None = None,
    en: bool = False,
    normalized: bool = False,
    E_norm: Sequence | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return differential Energy flux and its statistical error.

    Set ``normalized=True, E_norm=(10, 15)`` to return a spectrum relative to
    its mean 10--15 keV flux.
    """
    result = counts_nd(
        df,
        {"Energy": energy_bins},
        weights=weights,
        keepdims=("Energy",),
    )
    return counts_to_flux(
        result,
        exposure=exposure,
        geometric_factor=geometric_factor,
        en=en,
        normalized=normalized,
        E_norm=E_norm,
    )


def flux_t(
    df: pd.DataFrame,
    time_bins: Sequence,
    energy_range: Sequence,
    *,
    geometric_factor: float = 1.0,
    weights: ArrayLike | None = None,
    en: bool = False,
    normalized: bool = False,
    E_norm: Sequence | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return time-dependent flux in one Energy range.

    This function collapses Energy before returning, so Energy-band spectral
    normalization is unavailable. Use :func:`flux_t_e` instead.
    """
    if normalized or E_norm is not None:
        raise ValueError("flux_t has no retained Energy axis; use flux_t_e for normalized spectra.")
    if len(energy_range) != 2:
        raise ValueError("energy_range must contain exactly two edges.")
    time_edges = Time(time_bins).datetime
    exposure = np.diff(time_edges.asi8).astype(float) / 1.0e9
    if len(exposure) == 0 or np.any(~np.isfinite(exposure)) or np.any(exposure <= 0):
        raise ValueError("time_bins must define positive, increasing intervals.")
    result = counts_nd(
        df,
        {"Datetime": time_edges, "Energy": energy_range},
        weights=weights,
        keepdims=("Datetime",),
    )
    _, flux, flux_error = counts_to_flux(
        result,
        exposure=exposure,
        geometric_factor=geometric_factor,
        en=en,
    )
    return result.centers["Datetime"], flux, flux_error


def flux_t_e(
    df: pd.DataFrame,
    time_bins: Sequence,
    energy_bins: Sequence,
    *,
    geometric_factor: float = 1.0,
    weights: ArrayLike | None = None,
    en: bool = False,
    normalized: bool = False,
    E_norm: Sequence | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return time-Energy differential flux and its statistical error.

    With ``normalized=True``, every time-bin spectrum is separately divided
    by its mean differential flux within ``E_norm`` (keV).
    """
    time_edges = Time(time_bins).datetime
    exposure = np.diff(time_edges.asi8).astype(float) / 1.0e9
    if len(exposure) == 0 or np.any(~np.isfinite(exposure)) or np.any(exposure <= 0):
        raise ValueError("time_bins must define positive, increasing intervals.")
    result = counts_nd(
        df,
        {"Datetime": time_edges, "Energy": energy_bins},
        weights=weights,
        keepdims=("Datetime", "Energy"),
    )
    exposure_shape = [1] * result.raw_counts.ndim
    exposure_shape[result.dims.index("Datetime")] = len(exposure)
    energy, flux, flux_error = counts_to_flux(
        result,
        exposure=exposure.reshape(exposure_shape),
        geometric_factor=geometric_factor,
        en=en,
        normalized=normalized,
        E_norm=E_norm,
    )
    return result.centers["Datetime"], energy, flux, flux_error


def flux_to_temperature(
    combo_flux: tuple[ArrayLike, ArrayLike, ArrayLike],
    *,
    energy_axis: int = 0,
    min_energy: float = SCIME_MIN_ENERGY_KEV,
    max_energy: float | None = SCIME_MAX_ENERGY_KEV,
    sigma_func: Callable = sigma_HH,
    drift_velocity_kms: ArrayLike = 0.0,
    look_direction_cosine: ArrayLike = 1.0,
    particle_mass_kg: float = mass_H_kg,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fit N-dimensional differential flux and infer drift-aware temperature.

    For every spectrum along ``energy_axis``, fit the drifting-Maxwellian form
    ``ln(F(E) / (sigma(E) E)) = a - E_prime / T``, where
    ``E_prime = E + E_d - 2 sqrt(E E_d) mu sign(V_d)`` and
    ``E_d = m V_d**2 / 2``. ``drift_velocity_kms`` is the signed bulk speed
    along a chosen positive reference axis and ``look_direction_cosine`` is
    ``mu``, the cosine between the measured ENA velocity and that positive
    axis. Both may be scalars or arrays broadcastable to the output shape.
    Their defaults recover the original zero-drift Scime fit exactly.

    The Energy dimension is removed, so every returned array has one fewer
    dimension than flux. Return ``(temperature, temperature_error, slope,
    slope_error, r_squared)``. The returned slope is with respect to
    ``E_prime``. ``particle_mass_kg`` defaults to the proton mass.

    ``flux_error / flux`` is propagated to the logarithmic spectrum and used
    for inverse-variance weighted linear fitting. Temperature uncertainty uses
    first-order propagation through ``T = -1 / slope``. A non-negative slope
    leaves temperature and its error as ``NaN`` while retaining the fitted
    slope diagnostics. Spectra with fewer than three valid positive samples
    return ``NaN`` for every output. Drift velocity is treated as an externally
    known parameter; its uncertainty is not included in ``temperature_error``.
    A single ``mu`` per spectrum is exact for one look direction. An integrated
    FOV with a distribution of directions requires response-weighted forward
    modeling rather than replacing that distribution by its mean cosine.
    """
    if len(combo_flux) != 3:
        raise ValueError("combo_flux must contain (energy, flux, flux_error).")
    energy = np.asarray(combo_flux[0], dtype=float)
    flux = np.asarray(combo_flux[1], dtype=float)
    flux_err = np.asarray(combo_flux[2], dtype=float)

    if energy.ndim != 1:
        raise ValueError("energy must be one-dimensional.")
    if flux.ndim == 0:
        raise ValueError("flux must have at least one dimension.")
    if flux_err.shape != flux.shape:
        raise ValueError("flux_error must have the same shape as flux.")
    if not -flux.ndim <= energy_axis < flux.ndim:
        raise ValueError("energy_axis is out of bounds for flux.")
    energy_axis %= flux.ndim
    if flux.shape[energy_axis] != len(energy):
        raise ValueError("energy length must match flux along energy_axis.")
    particle_mass_value = np.asarray(particle_mass_kg, dtype=float)
    if particle_mass_value.ndim != 0 or not np.isfinite(particle_mass_value) or particle_mass_value <= 0:
        raise ValueError("particle_mass_kg must be a positive finite scalar.")
    particle_mass_kg = float(particle_mass_value)

    sigma = np.broadcast_to(np.asarray(sigma_func(energy), dtype=float), energy.shape)
    energy_valid = np.isfinite(energy) & np.isfinite(sigma)
    energy_valid &= (energy > 0) & (sigma > 0) & (energy >= min_energy)
    if max_energy is not None:
        energy_valid &= energy <= max_energy

    spectra = np.moveaxis(flux, energy_axis, 0)
    spectrum_errors = np.moveaxis(flux_err, energy_axis, 0)
    temperatures = np.full(spectra.shape[1:], np.nan, dtype=float)
    temperature_errors = np.full(spectra.shape[1:], np.nan, dtype=float)
    slopes = np.full(spectra.shape[1:], np.nan, dtype=float)
    slope_errors = np.full(spectra.shape[1:], np.nan, dtype=float)
    r_squared = np.full(spectra.shape[1:], np.nan, dtype=float)
    try:
        drift_velocities = np.broadcast_to(
            np.asarray(drift_velocity_kms, dtype=float),
            temperatures.shape,
        )
    except ValueError as exc:
        raise ValueError("drift_velocity_kms is not broadcastable to the temperature shape.") from exc
    try:
        direction_cosines = np.broadcast_to(
            np.asarray(look_direction_cosine, dtype=float),
            temperatures.shape,
        )
    except ValueError as exc:
        raise ValueError("look_direction_cosine is not broadcastable to the temperature shape.") from exc
    finite_cosines = direction_cosines[np.isfinite(direction_cosines)]
    if np.any(np.abs(finite_cosines) > 1.0 + 1e-12):
        raise ValueError("look_direction_cosine must lie between -1 and 1.")

    for index in np.ndindex(temperatures.shape):
        spectrum = spectra[(slice(None),) + index]
        spectrum_error = spectrum_errors[(slice(None),) + index]
        drift_velocity = drift_velocities[index]
        direction_cosine = direction_cosines[index]
        if not np.isfinite(drift_velocity):
            continue
        if drift_velocity == 0:
            direction_cosine = 0.0
        elif not np.isfinite(direction_cosine):
            continue
        valid = energy_valid & np.isfinite(spectrum) & (spectrum > 0)
        valid &= np.isfinite(spectrum_error) & (spectrum_error > 0)
        if np.count_nonzero(valid) < 3:
            continue
        measured_energy = energy[valid]
        drift_energy = (
            0.5 * particle_mass_kg * (abs(drift_velocity) * 1.0e3) ** 2
            / e_charge
            / 1.0e3
        )
        x = (
            measured_energy
            + drift_energy
            - 2.0
            * np.sqrt(measured_energy * drift_energy)
            * np.clip(direction_cosine, -1.0, 1.0)
            * np.sign(drift_velocity)
        )
        if np.unique(x).size < 2:
            continue
        y = np.log(spectrum[valid] / (sigma[valid] * measured_energy))
        y_error = spectrum_error[valid] / spectrum[valid]
        coefficients, covariance = np.polyfit(
            x,
            y,
            1,
            w=1.0 / y_error,
            cov="unscaled",
        )
        fitted_slope, fitted_intercept = coefficients
        fitted_slope_error = np.sqrt(covariance[0, 0])
        statistical_weights = 1.0 / y_error**2
        predicted = fitted_slope * x + fitted_intercept
        weighted_mean = np.average(y, weights=statistical_weights)
        ss_res = np.sum(statistical_weights * (y - predicted) ** 2)
        ss_tot = np.sum(statistical_weights * (y - weighted_mean) ** 2)

        slopes[index] = fitted_slope
        slope_errors[index] = fitted_slope_error
        if ss_tot > 0:
            r_squared[index] = 1.0 - ss_res / ss_tot
        if fitted_slope < 0:
            temperatures[index] = -1.0 / fitted_slope
            temperature_errors[index] = fitted_slope_error / fitted_slope**2
    return temperatures, temperature_errors, slopes, slope_errors, r_squared

"""Small, mission-independent helpers for relative energy spectra."""

import numpy as np


def normalize_flux_by_energy_band(
    energy,
    flux,
    flux_error=None,
    *,
    energy_axis=-1,
    energy_edges=None,
    E_norm,
):
    """Normalize each spectrum by its mean flux in an Energy reference band.

    ``energy`` contains the centres for the axis identified by
    ``energy_axis``.  When ``energy_edges`` is supplied, partially overlapping
    bins contribute in proportion to their overlap with ``E_norm``; otherwise
    all channel centres inside the closed reference interval have equal
    weights.  Every independent spectrum along the remaining axes receives
    its own reference value.

    If ``flux_error`` is supplied, its independent-count uncertainty is
    propagated through the ratio, including the covariance between a channel
    and the reference band that contains it.  The returned normalized flux is
    dimensionless and has a reference-band mean of one.
    """
    energy = np.asarray(energy, dtype=float)
    values = np.asarray(flux, dtype=float)
    if energy.ndim != 1:
        raise ValueError("energy must be a one-dimensional array of channel centres.")
    if values.ndim == 0:
        raise ValueError("flux must retain an Energy axis to be normalized.")
    if not isinstance(energy_axis, (int, np.integer)):
        raise TypeError("energy_axis must be an integer.")
    if not -values.ndim <= energy_axis < values.ndim:
        raise np.exceptions.AxisError(energy_axis, values.ndim)
    axis = energy_axis % values.ndim
    if values.shape[axis] != energy.size:
        raise ValueError("energy length must match flux along energy_axis.")

    try:
        lower, upper = np.asarray(E_norm, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("E_norm must contain exactly two finite Energy limits.") from exc
    if not np.isfinite([lower, upper]).all() or lower >= upper:
        raise ValueError("E_norm must contain two finite, increasing Energy limits.")

    if energy_edges is None:
        base_weights = ((energy >= lower) & (energy <= upper)).astype(float)
    else:
        edges = np.asarray(energy_edges, dtype=float)
        if edges.ndim != 1 or edges.size != energy.size + 1:
            raise ValueError("energy_edges must have exactly one more value than energy.")
        if not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0.0):
            raise ValueError("energy_edges must be finite and strictly increasing.")
        base_weights = np.clip(
            np.minimum(edges[1:], upper) - np.maximum(edges[:-1], lower),
            0.0,
            None,
        )
    if not np.any(base_weights > 0.0):
        raise ValueError("E_norm does not overlap any Energy channel.")

    moved_flux = np.moveaxis(values, axis, -1)
    shape = (1,) * (moved_flux.ndim - 1) + (energy.size,)
    weighted = base_weights.reshape(shape)
    valid = np.isfinite(moved_flux) & (weighted > 0.0)
    weight_sum = np.sum(np.where(valid, weighted, 0.0), axis=-1, keepdims=True)
    channel_weights = np.divide(
        weighted,
        weight_sum,
        out=np.zeros_like(moved_flux, dtype=float),
        where=valid & (weight_sum > 0.0),
    )
    reference = np.sum(
        np.where(valid, moved_flux * channel_weights, 0.0), axis=-1, keepdims=True
    )
    valid_reference = np.isfinite(reference) & (reference > 0.0)
    normalized_flux = np.divide(
        moved_flux,
        reference,
        out=np.full_like(moved_flux, np.nan, dtype=float),
        where=valid_reference,
    )

    if flux_error is None:
        return np.moveaxis(normalized_flux, -1, axis), None

    errors = np.asarray(flux_error, dtype=float)
    if errors.shape != values.shape:
        raise ValueError("flux_error must have the same shape as flux.")
    moved_error = np.moveaxis(errors, axis, -1)
    reference_variance = np.sum(
        (channel_weights * moved_error) ** 2, axis=-1, keepdims=True
    )
    covariance = channel_weights * moved_error**2
    normalized_variance = (
        moved_error**2 / reference**2
        + moved_flux**2 * reference_variance / reference**4
        - 2.0 * moved_flux * covariance / reference**3
    )
    normalized_error = np.sqrt(np.clip(normalized_variance, 0.0, None))
    normalized_error[~np.broadcast_to(valid_reference, normalized_error.shape)] = np.nan
    return (
        np.moveaxis(normalized_flux, -1, axis),
        np.moveaxis(normalized_error, -1, axis),
    )

"""Infer a time-independent relative ENA energy response from count data.

The response correction is shared by all supplied observations.  Each
observation may have one energy-independent normalization, but there is no
time-dependent or line-of-sight multiplier.  Consequently the inferred curve
can change spectral shape only; it cannot repair an incorrect time profile.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize, minimize_scalar


@dataclass(frozen=True)
class EnergyResponseEvent:
    """Observed and physical-model counts for one independent observation.

    All count arrays have shape ``(time, energy)``. ``source_counts`` excludes
    background and already includes the engineering detector response used by
    the physical forward model. ``background_counts`` is fixed independently.
    """

    name: str
    energy_keV: np.ndarray
    observed_counts: np.ndarray
    source_counts: np.ndarray
    background_counts: np.ndarray

    def __post_init__(self):
        energy = np.asarray(self.energy_keV, dtype=float)
        observed = np.asarray(self.observed_counts, dtype=float)
        source = np.asarray(self.source_counts, dtype=float)
        background = np.asarray(self.background_counts, dtype=float)
        if not self.name:
            raise ValueError("event name must be non-empty")
        if energy.ndim != 1 or energy.size < 3:
            raise ValueError("energy_keV must contain at least three energies")
        if np.any(~np.isfinite(energy)) or np.any(energy <= 0) or np.any(np.diff(energy) <= 0):
            raise ValueError("energy_keV must be finite, positive, and increasing")
        if observed.ndim != 2 or observed.shape[1] != energy.size:
            raise ValueError("observed_counts must have shape (time, energy)")
        if source.shape != observed.shape or background.shape != observed.shape:
            raise ValueError("source_counts and background_counts must match observed_counts")
        if any(np.any(~np.isfinite(item)) or np.any(item < 0) for item in (observed, source, background)):
            raise ValueError("count arrays must be finite and non-negative")
        if not np.any(source > 0):
            raise ValueError("source_counts must contain a positive value")
        object.__setattr__(self, "energy_keV", energy)
        object.__setattr__(self, "observed_counts", observed)
        object.__setattr__(self, "source_counts", source)
        object.__setattr__(self, "background_counts", background)


@dataclass(frozen=True)
class EnergyResponseFit:
    """Result of a shared relative-energy-response fit."""

    energy_keV: np.ndarray
    correction: np.ndarray
    anchor_energy_keV: float
    amplitudes: dict[str, float]
    predictions: dict[str, np.ndarray]
    poisson_deviance: float
    reduced_poisson_deviance: float
    objective: float
    success: bool
    message: str


def poisson_deviance(observed, expected):
    """Return the Poisson deviance, including zero-count cells."""

    observed = np.asarray(observed, dtype=float)
    expected = np.clip(np.asarray(expected, dtype=float), 1e-12, None)
    positive = observed > 0
    terms = expected - observed
    terms[positive] += observed[positive] * np.log(
        observed[positive] / expected[positive]
    )
    return float(2.0 * np.sum(terms))


def _anchor_weights(energy_keV, anchor_energy_keV):
    log_energy = np.log(np.asarray(energy_keV, dtype=float))
    log_anchor = np.log(float(anchor_energy_keV))
    if log_anchor < log_energy[0] or log_anchor > log_energy[-1]:
        raise ValueError("anchor_energy_keV must lie inside the energy grid")
    exact = np.flatnonzero(np.isclose(log_energy, log_anchor, rtol=0.0, atol=1e-12))
    weights = np.zeros(log_energy.size)
    if exact.size:
        weights[exact[0]] = 1.0
        return weights
    upper = int(np.searchsorted(log_energy, log_anchor))
    lower = upper - 1
    fraction = (log_anchor - log_energy[lower]) / (log_energy[upper] - log_energy[lower])
    weights[lower] = 1.0 - fraction
    weights[upper] = fraction
    return weights


def _curvature_matrix(energy_keV):
    x = np.log(np.asarray(energy_keV, dtype=float))
    matrix = np.zeros((x.size - 2, x.size))
    scale = np.median(np.diff(x)) ** 2
    for index in range(1, x.size - 1):
        left = x[index] - x[index - 1]
        right = x[index + 1] - x[index]
        factor = 2.0 / (left + right)
        matrix[index - 1, index - 1] = factor / left
        matrix[index - 1, index] = -factor * (1.0 / left + 1.0 / right)
        matrix[index - 1, index + 1] = factor / right
    return matrix * scale


def _validate_events(events):
    events = tuple(events)
    if not events:
        raise ValueError("at least one response-calibration event is required")
    names = [event.name for event in events]
    if len(set(names)) != len(names):
        raise ValueError("event names must be unique")
    reference = events[0].energy_keV
    for event in events[1:]:
        if event.energy_keV.shape != reference.shape or not np.allclose(
            event.energy_keV, reference, rtol=1e-10, atol=0.0
        ):
            raise ValueError("all events must use the same energy grid")
    return events


def apply_energy_response(event, correction, amplitude=1.0):
    """Apply a frozen response curve and scalar source normalization."""

    correction = np.asarray(correction, dtype=float)
    if correction.shape != event.energy_keV.shape:
        raise ValueError("correction must contain one value per energy")
    if np.any(~np.isfinite(correction)) or np.any(correction <= 0):
        raise ValueError("correction must be finite and positive")
    if not np.isfinite(amplitude) or amplitude <= 0:
        raise ValueError("amplitude must be finite and positive")
    return event.background_counts + amplitude * event.source_counts * correction[None, :]


def profile_event_amplitude(
    event,
    correction,
    *,
    amplitude_prior_sigma=np.log(3.0),
    log_amplitude_bounds=(-np.log(20.0), np.log(20.0)),
):
    """Fit only an energy-independent amplitude for a frozen response curve."""

    correction = np.asarray(correction, dtype=float)

    def objective(log_amplitude):
        prediction = apply_energy_response(event, correction, np.exp(log_amplitude))
        value = np.sum(
            prediction
            - event.observed_counts * np.log(np.clip(prediction, 1e-12, None))
        )
        if amplitude_prior_sigma is not None:
            value += 0.5 * (log_amplitude / amplitude_prior_sigma) ** 2
        return float(value)

    solution = minimize_scalar(
        objective,
        bounds=log_amplitude_bounds,
        method="bounded",
        options={"xatol": 1e-9},
    )
    return float(np.exp(solution.x))


def fit_shared_energy_response(
    events,
    *,
    anchor_energy_keV=10.0,
    smoothness_strength=4.0,
    response_prior_sigma=np.log(10.0),
    amplitude_prior_sigma=np.log(3.0),
    parameter_bounds=(-6.0, 6.0),
):
    """Fit one smooth ``F(E)`` shared by multiple observations.

    The log response is constrained to equal zero at ``anchor_energy_keV`` by
    log-linear interpolation between neighboring channels.  This removes the
    exact degeneracy between response normalization and per-event amplitudes.
    """

    events = _validate_events(events)
    energy = events[0].energy_keV
    anchor = _anchor_weights(energy, anchor_energy_keV)
    basis = null_space(anchor[None, :])
    curvature = _curvature_matrix(energy)
    n_response_parameters = basis.shape[1]
    initial = np.zeros(n_response_parameters + len(events))

    def unpack(parameters):
        log_response = basis @ parameters[:n_response_parameters]
        log_amplitudes = parameters[n_response_parameters:]
        return log_response, log_amplitudes

    def objective_and_gradient(parameters):
        log_response, log_amplitudes = unpack(parameters)
        response = np.exp(log_response)
        value = 0.0
        gradient_response = np.zeros(energy.size)
        gradient_amplitudes = np.zeros(len(events))
        for event_index, event in enumerate(events):
            signal = (
                np.exp(log_amplitudes[event_index])
                * event.source_counts
                * response[None, :]
            )
            prediction = np.clip(event.background_counts + signal, 1e-12, None)
            residual_factor = 1.0 - event.observed_counts / prediction
            value += np.sum(prediction - event.observed_counts * np.log(prediction))
            gradient_response += np.sum(residual_factor * signal, axis=0)
            gradient_amplitudes[event_index] = np.sum(residual_factor * signal)

        if smoothness_strength:
            curved = curvature @ log_response
            value += 0.5 * smoothness_strength * np.sum(curved**2)
            gradient_response += smoothness_strength * curvature.T @ curved
        if response_prior_sigma is not None:
            value += 0.5 * np.sum((log_response / response_prior_sigma) ** 2)
            gradient_response += log_response / response_prior_sigma**2
        if amplitude_prior_sigma is not None:
            value += 0.5 * np.sum((log_amplitudes / amplitude_prior_sigma) ** 2)
            gradient_amplitudes += log_amplitudes / amplitude_prior_sigma**2
        gradient = np.concatenate(
            (basis.T @ gradient_response, gradient_amplitudes)
        )
        return float(value), gradient

    solution = minimize(
        objective_and_gradient,
        initial,
        jac=True,
        method="L-BFGS-B",
        bounds=[parameter_bounds] * initial.size,
        options={"maxiter": 4000, "ftol": 1e-12, "gtol": 1e-7},
    )
    log_response, log_amplitudes = unpack(solution.x)
    correction = np.exp(log_response)
    amplitudes = {
        event.name: float(np.exp(log_amplitudes[index]))
        for index, event in enumerate(events)
    }
    predictions = {
        event.name: apply_energy_response(
            event, correction, amplitudes[event.name]
        )
        for event in events
    }
    deviance = sum(
        poisson_deviance(event.observed_counts, predictions[event.name])
        for event in events
    )
    cell_count = sum(event.observed_counts.size for event in events)
    fitted_count = n_response_parameters + len(events)
    degrees_of_freedom = max(cell_count - fitted_count, 1)
    return EnergyResponseFit(
        energy_keV=energy.copy(),
        correction=correction,
        anchor_energy_keV=float(anchor_energy_keV),
        amplitudes=amplitudes,
        predictions=predictions,
        poisson_deviance=float(deviance),
        reduced_poisson_deviance=float(deviance / degrees_of_freedom),
        objective=float(solution.fun),
        success=bool(solution.success),
        message=str(solution.message),
    )


def leave_one_event_out(events, **fit_kwargs):
    """Fit ``F(E)`` on all-but-one observations and score each held-out one.

    Both strict ``A=1`` and shape-only scores are returned.  The latter profiles
    one energy-independent amplitude on the held-out observation; it never fits
    an energy-dependent or time-dependent correction to held-out data.
    """

    events = _validate_events(events)
    if len(events) < 2:
        raise ValueError("leave-one-event-out validation requires at least two events")
    rows = []
    unity = np.ones_like(events[0].energy_keV)
    for held_out in events:
        training = [event for event in events if event.name != held_out.name]
        fit = fit_shared_energy_response(training, **fit_kwargs)
        baseline_amplitude = profile_event_amplitude(held_out, unity)
        corrected_amplitude = profile_event_amplitude(held_out, fit.correction)
        baseline_profiled = apply_energy_response(
            held_out, unity, baseline_amplitude
        )
        corrected_profiled = apply_energy_response(
            held_out, fit.correction, corrected_amplitude
        )
        baseline_strict = apply_energy_response(held_out, unity, 1.0)
        corrected_strict = apply_energy_response(held_out, fit.correction, 1.0)
        cell_count = held_out.observed_counts.size
        rows.append(
            {
                "held_out_event": held_out.name,
                "training_events": ",".join(event.name for event in training),
                "fit_success": fit.success,
                "baseline_profiled_amplitude": baseline_amplitude,
                "corrected_profiled_amplitude": corrected_amplitude,
                "baseline_profiled_reduced_deviance": poisson_deviance(
                    held_out.observed_counts, baseline_profiled
                ) / max(cell_count - 1, 1),
                "corrected_profiled_reduced_deviance": poisson_deviance(
                    held_out.observed_counts, corrected_profiled
                ) / max(cell_count - 1, 1),
                "baseline_strict_reduced_deviance": poisson_deviance(
                    held_out.observed_counts, baseline_strict
                ) / cell_count,
                "corrected_strict_reduced_deviance": poisson_deviance(
                    held_out.observed_counts, corrected_strict
                ) / cell_count,
                "trained_correction": fit.correction.copy(),
            }
        )
    return rows

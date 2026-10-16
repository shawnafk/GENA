import numpy as np

from libs.ena_response import (
    EnergyResponseEvent,
    apply_energy_response,
    fit_shared_energy_response,
    leave_one_event_out,
)


def synthetic_events():
    energy = np.array([5.0, 8.0, 12.0, 20.0, 40.0])
    log_response = np.log(np.array([0.72, 0.90, 1.12, 1.65, 3.0]))
    fraction = (np.log(10.0) - np.log(8.0)) / (np.log(12.0) - np.log(8.0))
    log_anchor = (1.0 - fraction) * log_response[1] + fraction * log_response[2]
    response = np.exp(log_response - log_anchor)
    base_spectrum = np.array([35.0, 46.0, 32.0, 15.0, 4.0])
    time_profiles = (
        np.array([0.4, 0.9, 1.3, 0.7]),
        np.array([1.2, 0.8, 0.5, 1.0]),
        np.array([0.6, 1.1, 0.9, 1.4]),
    )
    amplitudes = (0.75, 1.3, 1.8)
    events = []
    for index, (profile, amplitude) in enumerate(zip(time_profiles, amplitudes)):
        source = profile[:, None] * base_spectrum[None, :]
        background = np.full_like(source, 0.8 + 0.2 * index)
        observed = background + amplitude * source * response[None, :]
        events.append(
            EnergyResponseEvent(
                name=f"event-{index}",
                energy_keV=energy,
                observed_counts=observed,
                source_counts=source,
                background_counts=background,
            )
        )
    return events, response, amplitudes


def test_shared_response_recovers_synthetic_curve_and_event_amplitudes():
    events, response, amplitudes = synthetic_events()
    fit = fit_shared_energy_response(
        events,
        smoothness_strength=0.0,
        response_prior_sigma=None,
        amplitude_prior_sigma=None,
    )

    assert fit.success
    np.testing.assert_allclose(fit.correction, response, rtol=2e-5)
    np.testing.assert_allclose(
        [fit.amplitudes[event.name] for event in events], amplitudes, rtol=2e-5
    )
    assert fit.reduced_poisson_deviance < 1e-8


def test_response_anchor_is_log_interpolated_at_ten_kev():
    events, _, _ = synthetic_events()
    fit = fit_shared_energy_response(events)
    fraction = (np.log(10.0) - np.log(8.0)) / (np.log(12.0) - np.log(8.0))
    log_anchor = (
        (1.0 - fraction) * np.log(fit.correction[1])
        + fraction * np.log(fit.correction[2])
    )
    np.testing.assert_allclose(log_anchor, 0.0, atol=1e-12)


def test_leave_one_out_frozen_response_improves_held_out_spectral_shape():
    events, _, _ = synthetic_events()
    rows = leave_one_event_out(
        events,
        smoothness_strength=0.0,
        response_prior_sigma=None,
        amplitude_prior_sigma=None,
    )

    assert len(rows) == len(events)
    for row in rows:
        assert (
            row["corrected_profiled_reduced_deviance"]
            < row["baseline_profiled_reduced_deviance"]
        )


def test_response_changes_energy_shape_but_not_time_profile_within_a_channel():
    events, response, _ = synthetic_events()
    event = events[0]
    prediction = apply_energy_response(event, response, amplitude=2.0)
    source_after_background = prediction - event.background_counts
    for energy_index in range(event.energy_keV.size):
        np.testing.assert_allclose(
            source_after_background[:, energy_index]
            / event.source_counts[:, energy_index],
            2.0 * response[energy_index],
        )

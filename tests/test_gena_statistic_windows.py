from pathlib import Path
import runpy

import numpy as np
import pandas as pd


STATISTIC = runpy.run_path(
    Path(__file__).parents[1] / "gena" / "post" / "statistic.py"
)
hist = STATISTIC["hist"]
hist2d = STATISTIC["hist2d"]
window_ave = STATISTIC["window_ave"]
window_hist0d = STATISTIC["window_hist0d"]
window_hist1d = STATISTIC["window_hist1d"]
window_hist2d = STATISTIC["window_hist2d"]
flux_t = STATISTIC["flux_t"]
e_window = STATISTIC["get_bin_grid"]
flux_e = STATISTIC["flux_e"]
flux_t_e = STATISTIC["flux_t_e"]


def test_window_ave_groups_explicit_weights_with_their_events():
    events = pd.DataFrame(
        {
            "Datetime": pd.to_datetime(
                ["2025-01-01 00:00:00", "2025-01-01 00:00:30", "2025-01-01 00:02:00"]
            ),
            "Energy": [1.0, 2.0, 3.0],
        },
        index=[10, 11, 12],
    )
    weights = pd.Series([0.5, 1.5, 2.0], index=events.index)

    times, means, errors = window_ave(
        events,
        "Energy",
        window="1min",
        weights=weights,
    )

    pd.testing.assert_index_equal(
        pd.DatetimeIndex(times),
        pd.to_datetime(["2025-01-01 00:00:00", "2025-01-01 00:02:00"]),
    )
    np.testing.assert_allclose(means, [1.75, 3.0])
    np.testing.assert_allclose(errors[0], np.sqrt(0.5 / 1.6))
    assert np.isnan(errors[1])


def test_hist_returns_linear_and_log_bin_centers():
    linear_centers, linear_counts, _ = hist([0.5, 1.5], [0.0, 1.0, 2.0])
    log_centers, log_counts, _ = hist([1.5, 3.0], [1.0, 2.0, 4.0])

    np.testing.assert_allclose(linear_centers, [0.5, 1.5])
    np.testing.assert_array_equal(linear_counts, [1, 1])
    np.testing.assert_allclose(log_centers, [np.sqrt(2.0), np.sqrt(8.0)])
    np.testing.assert_array_equal(log_counts, [1, 1])


def test_hist2d_returns_both_axis_centers():
    x_centers, y_centers, counts, errors = hist2d(
        [1.5, 3.0],
        [0.5, 1.5],
        [[1.0, 2.0, 4.0], [0.0, 1.0, 2.0]],
    )

    np.testing.assert_allclose(x_centers, [np.sqrt(2.0), np.sqrt(8.0)])
    np.testing.assert_allclose(y_centers, [0.5, 1.5])
    np.testing.assert_array_equal(counts, [[1.0, 0.0], [0.0, 1.0]])
    np.testing.assert_array_equal(errors, counts)


def test_window_hist0d_returns_counts_and_keeps_empty_windows():
    events = pd.DataFrame(
        {
            "Datetime": pd.to_datetime(
                ["2025-01-01 00:00:00", "2025-01-01 00:00:30", "2025-01-01 00:02:00"]
            ),
        },
        index=[10, 11, 12],
    )
    weights = pd.Series([0.5, 1.5, 2.0], index=events.index)

    times, counts, errors = window_hist0d(events, window="1min")
    weighted_times, weighted_counts, weighted_errors = window_hist0d(
        events, window="1min", weights=weights
    )

    pd.testing.assert_index_equal(
        pd.DatetimeIndex(times),
        pd.to_datetime(
            ["2025-01-01 00:00:00", "2025-01-01 00:01:00", "2025-01-01 00:02:00"]
        ),
    )
    pd.testing.assert_index_equal(pd.DatetimeIndex(weighted_times), pd.DatetimeIndex(times))
    np.testing.assert_array_equal(counts, [2, 0, 1])
    np.testing.assert_allclose(errors, [np.sqrt(2), 0.0, 1.0])
    np.testing.assert_allclose(weighted_counts, [2.0, 0.0, 2.0])
    np.testing.assert_allclose(weighted_errors, [np.sqrt(2.5), 0.0, 2.0])


def test_window_hist1d_uses_group_indices_and_keeps_empty_windows():
    events = pd.DataFrame(
        {
            "Datetime": pd.to_datetime(
                ["2025-01-01 00:00:00", "2025-01-01 00:00:30", "2025-01-01 00:02:00"]
            ),
            "Energy": [0.5, 1.5, 2.5],
        },
        index=[10, 11, 12],
    )
    weights = pd.Series([0.5, 1.5, 2.0], index=events.index)

    times, centers, counts, errors = window_hist1d(
        events, "Energy", [0, 1, 2, 3], window="1min", weights=weights
    )

    assert len(times) == 3
    np.testing.assert_allclose(centers, [0.5, 1.5, 2.5])
    np.testing.assert_allclose(counts, [[0.5, 1.5, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 2.0]])
    np.testing.assert_allclose(errors, [[0.5, 1.5, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 2.0]])


def test_window_hist2d_uses_group_indices_and_keeps_empty_windows():
    events = pd.DataFrame(
        {
            "Datetime": pd.to_datetime(
                ["2025-01-01 00:00:00", "2025-01-01 00:00:30", "2025-01-01 00:02:00"]
            ),
            "Energy": [0.5, 0.6, 2.5],
            "Angle": [0.5, 0.6, 0.5],
        },
        index=[10, 11, 12],
    )
    weights = pd.Series([0.5, 1.5, 2.0], index=events.index)

    times, x_centers, y_centers, counts, errors = window_hist2d(
        events,
        ["Energy", "Angle"],
        [[0, 1, 2, 3], [0, 1, 2]],
        window="1min",
        weights=weights,
    )

    assert len(times) == 3
    np.testing.assert_allclose(x_centers, [0.5, 1.5, 2.5])
    np.testing.assert_allclose(y_centers, [0.5, 1.5])
    np.testing.assert_allclose(counts[0], [[2.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
    np.testing.assert_allclose(errors[0], [[np.sqrt(2.5), 0.0], [0.0, 0.0], [0.0, 0.0]])
    np.testing.assert_allclose(counts[1], np.zeros((3, 2)))
    np.testing.assert_allclose(errors[1], np.zeros((3, 2)))
    np.testing.assert_allclose(counts[2], [[0.0, 0.0], [0.0, 0.0], [2.0, 0.0]])
    np.testing.assert_allclose(errors[2], [[0.0, 0.0], [0.0, 0.0], [2.0, 0.0]])


def test_flux_normalizes_window_counts_and_errors_once():
    events = pd.DataFrame(
        {
            "Datetime": pd.to_datetime(
                [
                    "2025-01-01 00:00:00",
                    "2025-01-01 00:00:30",
                    "2025-01-01 00:01:00",
                    "2025-01-01 00:02:00",
                ]
            ),
            "Energy": [2.0, 4.0, 8.0, 3.0],
        },
        index=[10, 11, 12, 13],
    )
    weights = pd.Series([0.5, 1.5, 100.0, 2.0], index=events.index)
    denominator = (5.0 - 1.0) * 60.0 * 2.0 * 5.0
    time_bins = pd.to_datetime(
        [
            "2025-01-01 00:00:00",
            "2025-01-01 00:01:00",
            "2025-01-01 00:02:00",
            "2025-01-01 00:03:00",
        ]
    )

    times, particle_flux, particle_err = flux_t(
        events,
        time_bins,
        (1.0, 5.0),
        geometric_factor=10.0,
        weights=weights,
    )
    _, energy_flux, energy_err = flux_t(
        events,
        time_bins,
        (1.0, 5.0),
        en=True,
        geometric_factor=10.0,
        weights=weights,
    )

    assert len(times) == 3
    np.testing.assert_array_equal(
        times,
        np.array(
            [
                "2025-01-01T00:00:30",
                "2025-01-01T00:01:30",
                "2025-01-01T00:02:30",
            ],
            dtype="datetime64[ns]",
        ),
    )
    np.testing.assert_allclose(particle_flux, np.array([2.0, 0.0, 2.0]) / denominator)
    np.testing.assert_allclose(particle_err, np.array([np.sqrt(2.5), 0.0, 2.0]) / denominator)
    np.testing.assert_allclose(energy_flux, particle_flux * 3.0)
    np.testing.assert_allclose(energy_err, particle_err * 3.0)


def test_e_window_supports_one_energy_bin():
    centers, widths, bin_type = e_window([1.0, 3.0])

    np.testing.assert_allclose(centers, [2.0])
    np.testing.assert_allclose(widths, [2.0])
    assert bin_type == "linear"


def test_flux_e_normalizes_each_energy_bin_separately():
    events = pd.DataFrame({"Energy": [1.5, 3.0]})

    centers, values, errors = flux_e(
        events,
        [1.0, 2.0, 4.0],
        geometric_factor=10.0,
    )
    _, energy_values, energy_errors = flux_e(
        events,
        [1.0, 2.0, 4.0],
        geometric_factor=10.0,
        en=True,
    )
    _, exposure_values, exposure_errors = flux_e(
        events,
        [1.0, 2.0, 4.0],
        exposure=[2.0, 4.0],
        geometric_factor=10.0,
    )

    np.testing.assert_allclose(centers, [np.sqrt(2.0), np.sqrt(8.0)])
    np.testing.assert_allclose(values, [1.0 / 10.0, 1.0 / 20.0])
    np.testing.assert_allclose(errors, values)
    np.testing.assert_allclose(energy_values, values * centers)
    np.testing.assert_allclose(energy_errors, errors * centers)
    np.testing.assert_allclose(exposure_values, [1.0 / 20.0, 1.0 / 80.0])
    np.testing.assert_allclose(exposure_errors, exposure_values)


def test_flux_e_normalized_by_the_mean_reference_band_flux():
    events = pd.DataFrame({"Energy": [5.5] * 2 + [10.5] * 4 + [15.5] * 8})

    centers, values, errors = flux_e(
        events,
        [5.0, 10.0, 15.0, 20.0],
        normalized=True,
        E_norm=(10.0, 15.0),
    )

    np.testing.assert_allclose(centers, [7.5, 12.5, 17.5])
    np.testing.assert_allclose(values, [0.5, 1.0, 2.0])
    np.testing.assert_allclose(errors, [np.sqrt(0.1875), 0.0, np.sqrt(1.5)])


def test_flux_e_reference_band_uses_partial_energy_bin_overlap():
    events = pd.DataFrame({"Energy": [8.5] * 4 + [12.5] * 12})

    _, values, _ = flux_e(
        events,
        [8.0, 12.0, 16.0],
        normalized=True,
        E_norm=(10.0, 15.0),
    )

    reference_flux = (2.0 * 1.0 + 3.0 * 3.0) / 5.0
    np.testing.assert_allclose(values, [1.0 / reference_flux, 3.0 / reference_flux])


def test_flux_t_e_normalizes_every_time_spectrum_separately():
    events = pd.DataFrame(
        {
            "Datetime": pd.to_datetime([
                "2025-01-01 00:00:10",
                "2025-01-01 00:00:20",
                "2025-01-01 00:01:10",
                "2025-01-01 00:01:20",
                "2025-01-01 00:01:30",
            ]),
            "Energy": [5.5, 10.5, 5.5, 10.5, 10.5],
        }
    )
    bins = pd.to_datetime([
        "2025-01-01 00:00:00", "2025-01-01 00:01:00", "2025-01-01 00:02:00"
    ])

    _, _, values, errors = flux_t_e(
        events, bins, [5.0, 10.0, 15.0], normalized=True, E_norm=(10.0, 15.0)
    )

    np.testing.assert_allclose(values, [[1.0, 1.0], [0.5, 1.0]])
    np.testing.assert_allclose(errors, [[np.sqrt(2.0), 0.0], [np.sqrt(0.375), 0.0]])


def test_flux_t_e_returns_time_energy_flux_and_errors():
    events = pd.DataFrame(
        {
            "Datetime": pd.to_datetime(
                [
                    "2025-01-01 00:00:00",
                    "2025-01-01 00:00:30",
                    "2025-01-01 00:01:00",
                    "2025-01-01 00:02:00",
                ]
            ),
            "Energy": [1.5, 1.6, 3.0, 1.5],
        },
        index=[10, 11, 12, 13],
    )
    weights = pd.Series([0.5, 1.5, 2.0, 3.0], index=events.index)

    time_bins = pd.to_datetime(
        [
            "2025-01-01 00:00:00",
            "2025-01-01 00:01:00",
            "2025-01-01 00:02:00",
            "2025-01-01 00:03:00",
        ]
    )
    times, centers, values, errors = flux_t_e(
        events,
        time_bins,
        [1.0, 2.0, 4.0],
        geometric_factor=10.0,
        weights=weights,
    )
    _, _, energy_values, energy_errors = flux_t_e(
        events,
        time_bins,
        [1.0, 2.0, 4.0],
        geometric_factor=10.0,
        weights=weights,
        en=True,
    )

    assert len(times) == 3
    np.testing.assert_allclose(centers, [np.sqrt(2.0), np.sqrt(8.0)])
    np.testing.assert_allclose(values, [[2.0 / 600.0, 0.0], [0.0, 2.0 / 1200.0], [3.0 / 600.0, 0.0]])
    np.testing.assert_allclose(errors, [[np.sqrt(2.5) / 600.0, 0.0], [0.0, 2.0 / 1200.0], [3.0 / 600.0, 0.0]])
    np.testing.assert_allclose(energy_values, values * centers[None, :])
    np.testing.assert_allclose(energy_errors, errors * centers[None, :])


def test_counts_based_flux_retains_single_time_and_energy_bins():
    events = pd.DataFrame(
        {
            "Datetime": pd.to_datetime(["2025-01-01 00:00:10", "2025-01-01 00:00:20"]),
            "Energy": [1.5, 3.0],
        }
    )

    time_bins = pd.to_datetime(["2025-01-01 00:00:00", "2025-01-01 00:01:00"])
    centers, values, errors = flux_e(events, [1.0, 4.0])
    times, time_centers, time_values, time_errors = flux_t_e(
        events,
        time_bins,
        [1.0, 4.0],
    )

    np.testing.assert_allclose(centers, [2.5])
    np.testing.assert_allclose(values, [2.0 / 3.0])
    np.testing.assert_allclose(errors, [np.sqrt(2.0) / 3.0])
    assert len(times) == 1
    np.testing.assert_allclose(time_centers, centers)
    np.testing.assert_allclose(time_values, [[2.0 / 180.0]])
    np.testing.assert_allclose(time_errors, [[np.sqrt(2.0) / 180.0]])


def test_counts_based_flux_returns_all_explicit_empty_time_bins():
    events = pd.DataFrame(
        {
            "Datetime": pd.Series([], dtype="datetime64[ns]"),
            "Energy": pd.Series([], dtype=float),
        }
    )

    time_bins = pd.to_datetime(
        [
            "2025-01-01 00:00:00",
            "2025-01-01 00:01:00",
            "2025-01-01 00:03:00",
        ]
    )
    times, centers, values, errors = flux_t_e(
        events,
        time_bins,
        [1.0, 2.0, 4.0],
    )

    assert times.shape == (2,)
    np.testing.assert_allclose(centers, [np.sqrt(2.0), np.sqrt(8.0)])
    assert values.shape == errors.shape == (2, 2)
    np.testing.assert_array_equal(values, np.zeros((2, 2)))
    np.testing.assert_array_equal(errors, np.zeros((2, 2)))


def test_time_flux_uses_each_explicit_time_bin_width_as_exposure():
    events = pd.DataFrame(
        {
            "Datetime": pd.to_datetime(
                ["2025-01-01 00:00:30", "2025-01-01 00:02:00"]
            ),
            "Energy": [2.0, 2.0],
        }
    )
    time_bins = pd.to_datetime(
        [
            "2025-01-01 00:00:00",
            "2025-01-01 00:01:00",
            "2025-01-01 00:03:00",
        ]
    )

    _, values, errors = flux_t(events, time_bins, [1.0, 3.0])

    np.testing.assert_allclose(values, [1.0 / 120.0, 1.0 / 240.0])
    np.testing.assert_allclose(errors, values)


def test_time_flux_rejects_nonpositive_time_intervals():
    events = pd.DataFrame(
        {"Datetime": pd.to_datetime(["2025-01-01"]), "Energy": [2.0]}
    )
    repeated = pd.to_datetime(
        ["2025-01-01 00:00:00", "2025-01-01 00:01:00", "2025-01-01 00:01:00"]
    )

    with np.testing.assert_raises_regex(ValueError, "positive, increasing"):
        flux_t(events, repeated, [1.0, 3.0])
    with np.testing.assert_raises_regex(ValueError, "positive, increasing"):
        flux_t_e(events, repeated, [1.0, 2.0, 3.0])


def test_flux_t_rejects_energy_band_normalization_after_energy_is_collapsed():
    events = pd.DataFrame(
        {"Datetime": pd.to_datetime(["2025-01-01"]), "Energy": [2.0]}
    )
    with np.testing.assert_raises_regex(ValueError, "no retained Energy axis"):
        flux_t(
            events,
            pd.to_datetime(["2025-01-01 00:00:00", "2025-01-01 00:01:00"]),
            [1.0, 3.0],
            normalized=True,
            E_norm=(1.0, 2.0),
        )

import numpy as np
import pandas as pd
from gena.post.statistic import CountsNDResult, counts_nd


def test_counts_nd_uses_centers_matching_the_explicit_bin_spacing():
    frame = pd.DataFrame({"Energy": [0.0, 0.5, 1.5, 2.5, 3.0]})

    linear = counts_nd(frame, {"Energy": [0.0, 1.0, 2.0, 3.0]})
    logarithmic = counts_nd(pd.DataFrame({"Energy": [1.5, 3.0, 7.0]}), {"Energy": [1.0, 2.0, 4.0, 8.0]})
    irregular = counts_nd(pd.DataFrame({"Energy": [1.5, 3.0, 7.0]}), {"Energy": [1.0, 2.0, 5.0, 9.0]})

    np.testing.assert_array_equal(linear.raw_counts, [2, 1, 1])
    assert len(linear.centers["Energy"]) == len(linear.raw_counts) == 3
    np.testing.assert_allclose(linear.centers["Energy"], [0.5, 1.5, 2.5])
    np.testing.assert_allclose(logarithmic.centers["Energy"], np.sqrt([2.0, 8.0, 32.0]))
    np.testing.assert_allclose(irregular.centers["Energy"], [1.5, 3.5, 7.0])


def test_counts_nd_collapses_two_edge_axes_to_a_scalar():
    frame = pd.DataFrame(
        {
            "Energy": [1.0, 2.0, 3.0],
            "PX": [-1.0, 0.0, 1.0],
            "PY": [0.0, 0.0, 0.0],
        }
    )

    result = counts_nd(
        frame,
        {
            "Energy": [1.0, 3.0],
            "PX": [-2.0, 2.0],
            "PY": [-1.0, 1.0],
        },
    )

    assert result.dims == ()
    assert result.raw_counts.shape == ()
    assert result.raw_counts.item() == 2
    assert result.selected_size == 2
    assert set(result.collapsed_ranges) == {"Energy", "PX", "PY"}


def test_counts_nd_keepdims_retains_single_bin_axes():
    frame = pd.DataFrame(
        {
            "Datetime": pd.to_datetime(["2025-01-01T00:00:10", "2025-01-01T00:00:20"]),
            "Energy": [1.5, 3.0],
        }
    )
    result = counts_nd(
        frame,
        {
            "Datetime": pd.to_datetime(["2025-01-01T00:00:00", "2025-01-01T00:01:00"]),
            "Energy": [1.0, 4.0],
        },
        keepdims=("Datetime", "Energy"),
    )

    assert result.dims == ("Datetime", "Energy")
    assert result.raw_counts.shape == (1, 1)
    np.testing.assert_array_equal(result.raw_counts, [[2]])
    assert result.collapsed_ranges == {}


def test_counts_nd_rejects_misaligned_weights_and_unknown_keepdims():
    frame = pd.DataFrame({"Energy": [1.0, 2.0]})

    with np.testing.assert_raises_regex(ValueError, "aligned with df"):
        counts_nd(frame, {"Energy": [0.0, 3.0]}, weights=[1.0])
    with np.testing.assert_raises_regex(ValueError, "not present"):
        counts_nd(frame, {"Energy": [0.0, 3.0]}, keepdims=("Datetime",))


def test_counts_nd_preserves_requested_multidimensional_axis_order():
    frame = pd.DataFrame(
        {
            "Datetime": pd.to_datetime(
                [
                    "2025-01-01T00:00:00Z",
                    "2025-01-01T00:01:00Z",
                    "2025-01-01T00:02:00Z",
                    "2025-01-01T00:03:00Z",
                ]
            ),
            "Energy": [1.5, 2.5, 1.5, 2.5],
            "PX": [-0.5, 0.5, -0.5, 0.5],
            "PY": [0.0, 0.0, 0.0, 0.0],
            "X": [0.5, 0.5, 1.5, 1.5],
            "Y": [0.5, 1.5, 0.5, 1.5],
        }
    )
    time_edges = pd.to_datetime(
        [
            "2025-01-01T00:00:00Z",
            "2025-01-01T00:02:00Z",
            "2025-01-01T00:04:00Z",
        ]
    )

    result = counts_nd(
        frame,
        {
            "Datetime": time_edges,
            "Energy": [1.0, 3.0],
            "PX": [-1.0, 0.0, 1.0],
            "PY": [-1.0, 1.0],
            "X": [0.0, 1.0, 2.0],
            "Y": [0.0, 1.0, 2.0],
        },
    )

    assert isinstance(result, CountsNDResult)
    assert result.dims == ("Datetime", "PX", "X", "Y")
    assert result.raw_counts.shape == (2, 2, 2, 2)
    assert result.raw_counts.sum() == 4
    assert result.input_size == result.selected_size == 4
    np.testing.assert_array_equal(
        result.centers["Datetime"],
        np.array(["2025-01-01T00:01:00", "2025-01-01T00:03:00"], dtype="datetime64[ns]"),
    )


def test_counts_nd_normalizes_date_and_unix_timestamp_axes_identically():
    timestamps = np.array([1735689600.0, 1735689660.0, 1735689720.0])
    date_frame = pd.DataFrame({"Datetime": pd.to_datetime(timestamps, unit="s", utc=True)})
    timestamp_frame = pd.DataFrame({"Timestamp": timestamps})
    date_edges = ["2025-01-01T00:00:00Z", "2025-01-01T00:02:00Z", "2025-01-01T00:04:00Z"]
    timestamp_edges = [1735689600.0, 1735689720.0, 1735689840.0]

    by_date = counts_nd(date_frame, {"Datetime": date_edges})
    by_timestamp = counts_nd(timestamp_frame, {"Timestamp": timestamp_edges})

    np.testing.assert_array_equal(by_date.raw_counts, by_timestamp.raw_counts)
    np.testing.assert_allclose(by_timestamp.centers["Timestamp"], [1735689660.0, 1735689780.0])


def test_counts_nd_returns_raw_and_external_weighted_count_statistics():
    frame = pd.DataFrame({"X": [0.2, 0.4, 1.2]})
    weights = np.array([1.0, -2.0, 3.0])

    result = counts_nd(frame, {"X": [0.0, 1.0, 2.0]}, weights=weights)

    np.testing.assert_array_equal(result.raw_counts, [2, 1])
    np.testing.assert_allclose(result.raw_count_errors, [np.sqrt(2.0), 1.0])
    np.testing.assert_allclose(result.weighted_counts, [-1.0, 3.0])
    np.testing.assert_allclose(result.weighted_count_errors, [np.sqrt(5.0), 3.0])


def test_counts_nd_empty_selection_keeps_requested_shape():
    result = counts_nd(pd.DataFrame({"PX": [10.0]}), {"PX": [-1.0, 0.0, 1.0]})

    assert result.raw_counts.shape == (2,)
    np.testing.assert_array_equal(result.raw_counts, [0, 0])
    assert result.selected_size == 0

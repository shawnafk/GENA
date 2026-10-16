import datetime

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pytest
from PIL import Image

from plots import vis
from plots.basic import tplot


@pytest.mark.parametrize("arguments", [{}, {"f": None, "ax": None}])
def test_pass_axes_creates_axes_for_omitted_or_explicit_none(arguments):
    @vis.pass_axes()
    def plot(*, f=None, ax=None):
        ax.plot([0, 1])
        return f, ax

    figure, axis = plot(**arguments)
    try:
        assert axis.get_figure() is figure
        assert len(axis.lines) == 1
    finally:
        plt.close(figure)


def test_pass_axes_creates_grid_in_existing_figure():
    @vis.pass_axes((2, 2, True, False))
    def plot(f=None, ax=None):
        return f, ax

    figure = plt.figure()
    before = plt.get_fignums()
    try:
        result, axes = plot(f=figure, ax=None)
        assert result is figure
        assert axes == figure.axes
        assert len(axes) == 4
        assert axes[0].get_shared_x_axes().joined(axes[0], axes[3])
        assert not axes[0].get_shared_y_axes().joined(axes[0], axes[3])
        assert plt.get_fignums() == before
    finally:
        plt.close(figure)


@pytest.mark.parametrize("container", [lambda ax: ax, lambda ax: [ax],
                                       lambda ax: (ax,), lambda ax: np.array([[ax]])])
def test_pass_axes_reuses_single_axis_and_supports_positional_arguments(container):
    @vis.pass_axes()
    def plot(value, f=None, ax=None):
        ax.set_title(value)
        return f, ax

    figure, axis = plt.subplots()
    before = plt.get_fignums()
    try:
        result, reused = plot("example", None, container(axis))
        assert result is figure and reused is axis
        assert axis.get_title() == "example"
        assert plt.get_fignums() == before
    finally:
        plt.close(figure)


def test_pass_axes_validates_count_even_with_explicit_figure():
    @vis.pass_axes((2, 1, True, False))
    def plot(f=None, ax=None):
        return f, ax

    figure, axes = plt.subplots(2, 1)
    try:
        assert plot(ax=axes) == (figure, list(axes))
        with pytest.raises(ValueError, match="Expected 2 axes, got 1"):
            plot(f=figure, ax=axes[0])
        with pytest.raises(ValueError, match="Expected 2 axes, got 0"):
            plot(ax=[])
    finally:
        plt.close(figure)


def test_pass_axes_rejects_mixed_figures_and_invalid_axes():
    @vis.pass_axes((2, 1, False, False))
    def plot(f=None, ax=None):
        return f, ax

    first, axes = plt.subplots(2, 1)
    second, other = plt.subplots()
    try:
        with pytest.raises(ValueError, match="same supplied or inferred figure"):
            plot(ax=[axes[0], other])
        with pytest.raises(ValueError, match="same supplied or inferred figure"):
            plot(f=second, ax=axes)
        with pytest.raises(TypeError, match="Matplotlib Axes"):
            plot(ax=[1, 2])
    finally:
        plt.close(first)
        plt.close(second)


def test_savefig_creates_parent_and_uses_shared_defaults(tmp_path):
    figure, axis = plt.subplots()
    axis.plot([0, 1], [0, 1])
    path = tmp_path / "nested" / "figure.jpg"

    result = vis.savefig(figure, path)

    assert result == path
    with Image.open(path) as image:
        assert image.format == "JPEG"
        assert all(abs(value - 300) < 1 for value in image.info["dpi"])
    assert not plt.fignum_exists(figure.number)


def test_savefig_is_available_from_plots_namespace():
    from plots import savefig

    assert savefig is vis.savefig


def test_tplot_is_the_single_entrypoint_for_dataset_and_xy_data():
    dataset = {
        "density": {
            "x": ["2025-01-01T00:00:00", "2025-01-01T00:01:00"],
            "y": [1.0, 2.0],
        }
    }
    figure, axes = plt.subplots(2, 1)

    tplot(dataset, key="density", ax=axes[0])
    tplot(dataset["density"], ax=axes[1])

    assert axes[0].lines[0].get_label() == "density"
    np.testing.assert_allclose(axes[0].lines[0].get_ydata(), [1.0, 2.0])
    np.testing.assert_allclose(axes[1].lines[0].get_ydata(), [1.0, 2.0])
    plt.close(figure)


def test_sat_coord_formatter_only_formats_precomputed_position():
    tick = mdates.date2num(
        datetime.datetime(2025, 1, 2, 3, 4, tzinfo=datetime.timezone.utc)
    )

    label = vis._sat_coord_formatter(tick, np.array([1.25, -2.5, 3.75]))

    assert label == "02 03:04\n1.2\n-2.5\n3.8"


def test_format_pos_date_caches_ephemeris_lookup_for_each_tick(monkeypatch):
    calls = []

    def fake_satellites(timestamps, sat_name):
        timestamps = np.asarray(timestamps)
        calls.append((timestamps.copy(), sat_name))
        return np.column_stack(
            (
                np.arange(timestamps.size),
                np.arange(timestamps.size) + 10,
                np.arange(timestamps.size) + 20,
            )
        )

    monkeypatch.setattr(vis.epms, "satellites", fake_satellites)
    figure, axis = plt.subplots()
    axis.set_xlim(
        datetime.datetime(2025, 1, 1, tzinfo=datetime.timezone.utc),
        datetime.datetime(2025, 1, 2, tzinfo=datetime.timezone.utc),
    )

    vis.format_pos_date(axis, sat_name="MMS1")
    formatter = axis.xaxis.get_major_formatter()
    locations = axis.xaxis.get_majorticklocs()
    labels = formatter.format_ticks(locations)

    assert locations.size <= 8
    assert len(calls) == locations.size
    assert all(call[1] == "MMS1" for call in calls)
    assert all(call[0].shape == () for call in calls)
    assert "0.0\n10.0\n20.0" in labels[0]

    headings = [
        artist for artist in axis.texts
        if artist.get_gid() == "satellite-position-row-heading"
    ]
    assert len(headings) == 1
    assert headings[0].get_text() == "T\nX\nY\nZ"

    figure.canvas.draw()
    visible_labels = [label for label in axis.get_xticklabels() if label.get_visible()]
    assert visible_labels[0].get_ha() == "left"
    assert visible_labels[-1].get_ha() == "right"
    heading_anchor = headings[0].get_transform().transform(headings[0].get_position())
    tick_anchor = visible_labels[0].get_transform().transform(visible_labels[0].get_position())
    assert heading_anchor[1] == tick_anchor[1]

    formatter.format_ticks(locations)
    assert len(calls) == locations.size

    formatter(locations[-1] + 1.0)
    assert len(calls) == locations.size + 1
    plt.close(figure)

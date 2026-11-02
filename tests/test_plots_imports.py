"""Import and public-API guardrails for the flattened plotting package."""

from __future__ import annotations

import importlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_public_api_imports():
    import plots

    from plots import Panel, env, gena, joint, mms, omni, overview, plot_region_code
    from plots import labels, panels, savefig, scene, themis, tplot

    assert callable(tplot)
    assert callable(savefig)
    assert callable(overview)
    assert callable(env)
    assert callable(plot_region_code)
    assert Panel is joint.Panel
    assert callable(labels.describe)
    assert callable(panels.plot_panel)
    assert callable(scene.prepare)
    for module in (gena, themis, mms, omni):
        assert module.__name__.startswith("plots.")


def test_lightweight_modules_import():
    import common.frame  # noqa: F401
    import gena.post.export  # noqa: F401
    import gena.utils  # noqa: F401
    import plots.ani  # noqa: F401
    import plots.v3d  # noqa: F401


def test_removed_legacy_paths_are_gone():
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("plots.vis")
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("plots.basic")
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("plots.advance")


def test_import_does_not_mutate_rcparams():
    """Importing plotting modules must not apply a global style implicitly."""
    script = (
        "import matplotlib\n"
        "before = dict(matplotlib.rcParams)\n"
        "import plots.gena, plots.themis, plots.omni\n"
        "after = dict(matplotlib.rcParams)\n"
        "assert before == after, 'importing plotting modules changed rcParams'\n"
    )
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(REPO_ROOT)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="gena-mpl-"))
    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr or completed.stdout


def test_theme_is_explicit_and_reversible():
    import plots.theme

    original = dict(matplotlib.rcParams)
    try:
        updated = plots.theme.apply_theme({"font.size": 17})
        assert updated["font.size"] == 17
        assert matplotlib.rcParams["font.size"] == 17
    finally:
        matplotlib.rcParams.update(original)


def test_import_does_not_load_geopack():
    """3-D imports must not trigger geopack's IGRF network update."""
    script = (
        "import sys\n"
        "import plots.v3d\n"
        "import libs.space\n"
        "assert 'geopack.geopack' not in sys.modules, 'geopack was imported eagerly'\n"
    )
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(REPO_ROOT)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="gena-mpl-"))
    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr or completed.stdout


def test_tplot_accepts_plain_dataframe():
    from plots import tplot

    frame = pd.DataFrame(
        {
            "Datetime": pd.date_range("2025-12-04T00:00:00", periods=4, freq="min"),
            "value": [1.0, 2.0, 3.0, 4.0],
        }
    )
    figure, axis = tplot(frame, "value")
    try:
        assert len(axis.lines) == 1
    finally:
        plt.close(figure)


def test_variable_labels_feed_tplot_axis():
    from plots import tplot
    from plots.labels import axis_label, describe

    info = describe("thb_peir_density")
    assert info.unit == "cm$^{-3}$"
    assert axis_label(info) == r"$n_i$ (cm$^{-3}$)"

    frame = pd.DataFrame(
        {
            "Datetime": pd.date_range("2025-12-04T00:00:00", periods=3, freq="min"),
            "BZ_GSM": [1.0, 2.0, 3.0],
        }
    )
    figure, axis = tplot(frame, "BZ_GSM")
    try:
        assert axis.get_ylabel() == r"$B_z$ (nT)"
        assert axis.lines[0].get_label() == r"$B_z$"
    finally:
        plt.close(figure)


def test_mms_overview_accepts_one_product():
    from plots import mms

    times = pd.date_range("2025-12-04T00:00:00", periods=4, freq="min")
    frame = pd.DataFrame({"Datetime": times})
    for name in (
        "numberdensity_fast",
        "tempperp_fast",
        "temppara_fast",
        "bulkv_gse_fast",
    ):
        frame[name] = [1.0, 2.0, 3.0, 4.0]

    figure, axes = mms.overview(fpi_data=frame)
    try:
        assert len(axes) == 3
    finally:
        plt.close(figure)


def test_themis_overview_accepts_partial_products():
    from plots import themis

    times = pd.date_range("2025-12-04T00:00:00", periods=4, freq="min")
    fgm = {
        "thb_fgs_btotal": {"x": times, "y": [1.0, 2.0, 3.0, 4.0]},
        "thb_fgs_gsm": {"x": times, "y": [[1.0, 0.0, 0.0]] * 4},
    }
    figure, axes = themis.overview(fgm=fgm)
    try:
        assert len(axes) == 1
    finally:
        plt.close(figure)


def test_panel_spec_scalar():
    from plots.panels import PanelSpec, plot_panel

    times = pd.date_range("2025-12-04T00:00:00", periods=4, freq="min")
    frame = pd.DataFrame({"Datetime": times, "BZ_GSM": [1.0, 2.0, 3.0, 4.0]})
    figure, axis = plot_panel(frame, PanelSpec(("BZ_GSM",), zero=True))
    try:
        assert axis.get_ylabel() == r"$B_z$ (nT)"
        assert len(axis.lines) == 2  # data + zero line
    finally:
        plt.close(figure)


def test_themis_overview_with_esa_spectrum():
    import numpy as np

    from plots import themis

    times = pd.date_range("2025-12-04T00:00:00", periods=4, freq="min")
    esa = {
        "thb_peir_density": {"x": times, "y": [1.0, 2.0, 3.0, 4.0]},
        "thb_peir_avgtemp": {"x": times, "y": [10.0, 20.0, 30.0, 40.0]},
        "thb_peir_en_eflux": {
            "x": times,
            "y": np.ones((4, 3)),
            "v": np.array([100.0, 1000.0, 10000.0]),
        },
    }
    figure, axes = themis.overview(esa=esa)
    try:
        assert len(axes) == 3
    finally:
        plt.close(figure)


def test_plasma_mach_beta_helpers():
    import numpy as np

    from libs import plasma

    magnetic = np.array([1.0, 10.0])
    density = np.array([1.0, 5.0])
    velocity = np.array([100.0, 300.0])
    temperature = np.array([100.0, 500.0])
    assert np.allclose(
        plasma.alfven_mach(magnetic, density, velocity),
        velocity * 1000.0 / plasma.alfven_velocity(magnetic * 1e-5, density, 1),
    )
    assert np.allclose(
        plasma.plasma_beta(magnetic, density, temperature),
        plasma.beta(density, temperature, magnetic * 1e-5),
    )

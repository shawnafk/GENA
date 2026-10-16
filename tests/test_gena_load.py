import importlib

import numpy as np
import pandas as pd
import pytest

from gena.frame import ParameterDataFrame, ParticleDataFrame


gena_load = importlib.import_module("load.load")


def _write(path, rows):
    pd.DataFrame(rows).to_parquet(path, engine="pyarrow", index=False)


def _configure_pqt(monkeypatch, path):
    monkeypatch.setattr(
        gena_load,
        "get_configured_paths",
        lambda: {"pqt_base": str(path)},
    )


def test_single_pdf_file_is_sorted_deduplicated_and_typed(tmp_path):
    path = tmp_path / "20240101_pdf.pqt"
    _write(
        path,
        {
            "Timestamp": [2.0, 1.0, 1.0],
            "Value": [20, 10, 99],
        },
    )

    result = gena_load.load_gena(file=path)

    assert isinstance(result, ParameterDataFrame)
    assert result["Timestamp"].tolist() == [1.0, 2.0]
    assert result["Value"].tolist() == [10, 20]
    assert result.index.tolist() == [0, 1]
    assert result["Datetime"].dtype == np.dtype("datetime64[ns]")


def test_file_and_range_sources_are_mutually_exclusive(tmp_path):
    path = tmp_path / "20240101_pdf.pqt"
    _write(path, {"Timestamp": [1.0, 2.0, 3.0, 4.0], "Value": range(4)})

    with pytest.raises(ValueError, match="exactly one"):
        gena_load.load_gena()
    with pytest.raises(ValueError, match="exactly one"):
        gena_load.load_gena(file=path, trange=[1.0, 2.0], product="pdf")
    with pytest.raises(ValueError, match="product cannot"):
        gena_load.load_gena(file=path, product="pdf")
    with pytest.raises(ValueError, match="product is required"):
        gena_load.load_gena(trange=[1.0, 2.0])
    with pytest.raises(TypeError):
        gena_load.load_gena([1.0, 2.0], "pdf")


def test_range_mode_discovers_nested_daily_files_once(tmp_path, monkeypatch):
    _configure_pqt(monkeypatch, tmp_path)
    nested = tmp_path / "nested"
    nested.mkdir()
    _write(tmp_path / "20240101_pdf.pqt", {"Timestamp": [1704067200.0], "Value": [1]})
    _write(nested / "20240102_pdf.pqt", {"Timestamp": [1704153600.0], "Value": [2]})
    _write(nested / "20240103_pdf.pqt", {"Timestamp": [1704240000.0], "Value": [3]})

    result = gena_load.load_gena(
        trange=["2024-01-01T00:00:00Z", "2024-01-02T00:00:00Z"],
        product="pdf",
    )

    assert result["Value"].tolist() == [1, 2]


def test_particle_products_keep_discrete_pose_columns_without_pose_methods(tmp_path):
    path = tmp_path / "20240101_l2_a.pqt"
    _write(
        path,
        {
            "UUID": ["a", "b", "c"],
            "Timestamp": [1.0, 1.0, 1.0],
            "Value": [10, 10, 20],
            "MJ2000_Star_X": [1.0, 2.0, 3.0],
            "MJ2000_q1": [0.1, 0.2, 0.3],
        },
    )

    result = gena_load.load_gena(file=path, columns=["Timestamp", "MJ2000_Star_X", "MJ2000_q1"])

    assert isinstance(result, ParticleDataFrame)
    assert result["Timestamp"].tolist() == [1.0, 1.0, 1.0]
    assert result["MJ2000_Star_X"].tolist() == [1.0, 2.0, 3.0]
    assert not hasattr(result, "pos")
    assert not hasattr(result, "rot")
    assert not hasattr(result, "pos_at")
    assert not hasattr(result, "rot_at")


def test_empty_range_result_uses_configured_path(tmp_path, monkeypatch):
    _configure_pqt(monkeypatch, tmp_path)
    result = gena_load.load_gena(
        trange=["2024-01-01", "2024-01-01 23:59:59"],
        product="pdf",
    )
    assert isinstance(result, ParameterDataFrame)
    assert result.empty


def _parameter_frame():
    return ParameterDataFrame(
        {
            "Datetime": pd.to_datetime([0.0, 1.0, 2.0], unit="s"),
            "Timestamp": [0.0, 1.0, 2.0],
            "MJ2000_Star_X": [1.0, 2.0, 3.0],
            "MJ2000_Star_Y": [2.0, 3.0, 4.0],
            "MJ2000_Star_Z": [3.0, 4.0, 5.0],
            "MJ2000_q1": [0.1, 0.15, 0.2],
            "MJ2000_q2": [0.1, 0.1, 0.1],
            "MJ2000_q3": [0.1, 0.1, 0.1],
        }
    )


def test_parameter_position_interpolation_preserves_input_order():
    result = _parameter_frame().pos_at([2.0, 0.0, 1.0])

    np.testing.assert_allclose(
        result,
        [[3.0, 4.0, 5.0], [1.0, 2.0, 3.0], [2.0, 3.0, 4.0]],
    )


def test_parameter_rotation_interpolation_preserves_input_order():
    pdf = _parameter_frame()
    _, source_rot = pdf._rot_from_pdf()

    result = pdf.rot_at([2.0, 0.0, 1.0])

    np.testing.assert_allclose(
        result.as_quat(),
        source_rot.as_quat()[[2, 0, 1]],
        atol=1e-12,
    )


def test_parameter_interpolation_has_explicit_empty_and_coverage_behavior():
    pdf = _parameter_frame()

    assert pdf.pos_at([]).shape == (0, 3)
    with pytest.raises(ValueError, match="at least one"):
        pdf.rot_at([])
    with pytest.raises(ValueError, match="outside PDF coverage"):
        pdf.pos_at([3.0])
    with pytest.raises(ValueError, match="outside PDF coverage"):
        pdf.rot_at([3.0])


@pytest.mark.parametrize("gap", [0, 1.5])
def test_pdf_pose_cleaners_require_a_positive_integer_gap(gap):
    pdf = _parameter_frame()

    with pytest.raises(ValueError, match="positive integer"):
        pdf._loc_from_pdf(gap=gap)
    with pytest.raises(ValueError, match="positive integer"):
        pdf._rot_from_pdf(gap=gap)

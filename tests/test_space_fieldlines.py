import numpy as np
import pytest

from libs import space


def test_t96_y0_field_lines_preserve_xyz_and_model_parameters(monkeypatch):
    recalculated = []
    transformed = []
    traced = []

    def fake_recalc(ut):
        recalculated.append(ut)

    def fake_trace(x, y, z, direction, **kwargs):
        traced.append((x, y, z, direction, kwargs))
        points = np.array(
            [
                [x, x + direction],
                [y, y + 0.25 * direction],
                [z, z + direction],
            ]
        )
        return x, y, z, *points

    def fake_smgsm(x, y, z, direction):
        transformed.append((x, y, z, direction))
        return x, y, z

    monkeypatch.setattr(space.geopack, "recalc", fake_recalc)
    monkeypatch.setattr(space.geopack, "smgsm", fake_smgsm)
    monkeypatch.setattr(space.geopack, "trace", fake_trace)
    model = space.T96_fieldlines(
        1_742_000_000.0,
        2.5,
        -35.0,
        4.0,
        -6.0,
        Nlat=3,
    )

    lines = model.trace_lines_y0(
        seed_radius=1.05,
        footprint_latitude_range=(70.0, 80.0),
        rlim=80,
    )

    assert recalculated == [1_742_000_000.0]
    assert len(lines) == 3
    assert all(line.shape == (3, 2) for line in lines)
    assert len(transformed) == 3
    assert all(call[0] < 0 and call[1] == 0 and call[2] > 0 for call in transformed)
    assert all(call[3] == 1 for call in transformed)
    assert all(abs(call[1]) < 1e-12 for call in traced)
    assert all(call[3] == 1 for call in traced)
    assert all(call[4]["exname"] == "t96" for call in traced)
    assert all(call[4]["rlim"] == 80 for call in traced)
    assert all(call[4]["r0"] == 1 for call in traced)
    assert all(call[4]["maxloop"] == 2000 for call in traced)
    for call in traced:
        np.testing.assert_allclose(call[4]["parmod"], [2.5, -35.0, 4.0, -6.0])


def test_t96_y0_field_lines_validate_model_inputs():
    model = space.T96_fieldlines(0.0, np.nan, -20.0, 0.0, -5.0)

    try:
        model.trace_lines_y0()
    except ValueError as exc:
        assert "finite" in str(exc)
    else:
        raise AssertionError("non-finite T96 parameters should fail")


def test_t96_neutral_sheet_finds_total_bx_zero_surface(monkeypatch):
    recalculated = []
    external_calls = []

    def fake_recalc(ut):
        recalculated.append(ut)
        return 0.2

    def fake_external(name, parameters, psi, x, y, z):
        external_calls.append((name, parameters.copy(), psi, x, y, z))
        return z - 0.1 * x + 0.2 * y, 0.0, 0.0

    monkeypatch.setattr(space.geopack, "recalc", fake_recalc)
    monkeypatch.setattr(space.geopack, "call_external_model", fake_external)
    monkeypatch.setattr(
        space.geopack,
        "igrf_gsm",
        lambda x, y, z: (0.0, 0.0, 0.0),
    )
    sheet = space.T96NeutralSheet(
        1_742_000_000.0,
        2.5,
        -35.0,
        4.0,
        -6.0,
        x_range=(-10.0, -6.0),
        y_range=(-2.0, 2.0),
        z_range=(-5.0, 5.0),
        nx=3,
        ny=3,
        nz=7,
        half_thickness=2.0,
    )

    X, Y, Z = sheet.surface_xyz()
    lower, upper = sheet.thickness_xyz()

    np.testing.assert_allclose(Z, 0.1 * X - 0.2 * Y)
    np.testing.assert_allclose(lower[2], Z - 2.0)
    np.testing.assert_allclose(upper[2], Z + 2.0)
    assert sheet.z_at(-8.0, 0.0) == pytest.approx(-0.8)
    assert np.isnan(sheet.z_at(-12.0, 0.0))
    assert recalculated == [1_742_000_000.0]
    assert external_calls
    assert all(call[0] == "t96" for call in external_calls)
    for call in external_calls:
        np.testing.assert_allclose(call[1], [2.5, -35.0, 4.0, -6.0])

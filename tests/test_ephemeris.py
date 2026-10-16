import numpy as np
import pytest
from xml.etree import ElementTree as ET

from libs import ephemeris


def test_satellites_reuses_cached_interpolators(monkeypatch):
    load_calls = []
    interpolation_calls = []
    orbit = np.array(
        [
            [0.0, 0.0, 10.0, 20.0],
            [1.0, 1.0, 11.0, 21.0],
            [2.0, 2.0, 12.0, 22.0],
            [3.0, 3.0, 13.0, 23.0],
        ]
    )

    class FakeOrbitFile:
        def __enter__(self):
            return {"MMS1": orbit}

        def __exit__(self, exc_type, exc_value, traceback):
            return False

    def fake_load(*args, **kwargs):
        load_calls.append((args, kwargs))
        return FakeOrbitFile()

    def fake_interp1d(x, y, **kwargs):
        interpolation_calls.append((x.copy(), y.copy(), kwargs))
        return lambda values: np.interp(values, x, y)

    monkeypatch.setattr(ephemeris.np, "load", fake_load)
    monkeypatch.setattr(ephemeris, "interp1d", fake_interp1d)
    ephemeris._satellite_interpolators.cache_clear()
    try:
        first = ephemeris.satellites([0.25], "MMS1")
        second = ephemeris.satellites([0.75], "MMS1")
    finally:
        ephemeris._satellite_interpolators.cache_clear()

    assert len(load_calls) == 1
    assert len(interpolation_calls) == 3
    np.testing.assert_allclose(first, [[0.25, 10.25, 20.25]])
    np.testing.assert_allclose(second, [[0.75, 10.75, 20.75]])


def test_satellites_rejects_extrapolation_by_default(monkeypatch):
    orbit = np.array(
        [
            [0.0, 0.0, 10.0, 20.0],
            [1.0, 1.0, 11.0, 21.0],
            [2.0, 2.0, 12.0, 22.0],
            [3.0, 3.0, 13.0, 23.0],
        ]
    )

    class FakeOrbitFile:
        def __enter__(self):
            return {"MMS1": orbit}

        def __exit__(self, exc_type, exc_value, traceback):
            return False

    monkeypatch.setattr(ephemeris.np, "load", lambda *args, **kwargs: FakeOrbitFile())
    ephemeris._satellite_interpolators.cache_clear()
    try:
        with pytest.raises(ValueError, match="ephemeris covers"):
            ephemeris.satellites([4.0], "MMS1")
        extrapolated = ephemeris.satellites([4.0], "MMS1", extrapolate=True)
    finally:
        ephemeris._satellite_interpolators.cache_clear()

    assert extrapolated.shape == (1, 3)


def test_download_ssc_locations_converts_km_to_earth_radii(monkeypatch):
    xml = b"""<Response xmlns="http://sscweb.gsfc.nasa.gov/schema"><Result><StatusCode>Success</StatusCode><StatusSubCode>Success</StatusSubCode><Data><Id>mms1</Id><Coordinates><CoordinateSystem>Gse</CoordinateSystem><X>6378.1</X><Y>0</Y><Z>0</Z></Coordinates><Time>2026-04-07T00:00:00Z</Time></Data></Result></Response>"""
    monkeypatch.setattr(ephemeris, "_ssc_xml", lambda *args, **kwargs: ET.fromstring(xml))

    result = ephemeris._download_ssc_locations(
        "mms1",
        "2026-04-07T00:00:00Z",
        "2026-04-07T00:00:00Z",
        resolution_factor=16,
    )

    np.testing.assert_allclose(result[:, 1:], [[1.0, 0.0, 0.0]])


def test_update_satellite_orbits_merges_and_writes_metadata(tmp_path, monkeypatch):
    orbit_path = tmp_path / "orbits.npz"
    metadata_path = tmp_path / "orbits.meta.json"
    np.savez(
        orbit_path,
        MMS1=np.array(
            [
                [0.0, 1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0, 0.0],
                [2.0, 1.0, 0.0, 0.0],
                [3.0, 1.0, 0.0, 0.0],
            ]
        ),
    )
    monkeypatch.setattr(
        ephemeris,
        "_ssc_observatory_coverage",
        lambda **kwargs: {
            "mms1": (
                ephemeris.pd.Timestamp(0, unit="s", tz="UTC"),
                ephemeris.pd.Timestamp(10, unit="s", tz="UTC"),
            )
        },
    )
    monkeypatch.setattr(
        ephemeris,
        "_download_ssc_interval",
        lambda *args, **kwargs: np.array(
            [
                [4.0, 2.0, 0.0, 0.0],
                [5.0, 3.0, 0.0, 0.0],
            ]
        ),
    )

    coverage = ephemeris.update_satellite_orbits(
        "MMS1",
        end=ephemeris.pd.Timestamp(5, unit="s", tz="UTC"),
        orbit_path=orbit_path,
        metadata_path=metadata_path,
    )

    with np.load(orbit_path) as archive:
        updated = archive["MMS1"]
    assert updated[-1, 0] == 5.0
    assert updated[-1, 1] == 3.0
    assert coverage["MMS1"][1] == ephemeris.pd.Timestamp(5, unit="s", tz="UTC")
    assert '"coordinate_system": "GSE"' in metadata_path.read_text()


def test_satellites_update_extends_cache_before_interpolation(monkeypatch):
    class FakeInterpolators:
        updated = False

        def __call__(self, key):
            stop = 3.0 if self.updated else 1.0
            return (
                lambda values: np.asarray(values, dtype=float),
                lambda values: np.zeros_like(np.asarray(values, dtype=float)),
                lambda values: np.zeros_like(np.asarray(values, dtype=float)),
                0.0,
                stop,
            )

        def cache_clear(self):
            return None

    interpolators = FakeInterpolators()

    def fake_update(key, *, end):
        interpolators.updated = True

    monkeypatch.setattr(ephemeris, "_satellite_interpolators", interpolators)
    monkeypatch.setattr(ephemeris, "update_satellite_orbits", fake_update)

    result = ephemeris.satellites([2.0], "MMS1", update=True)

    np.testing.assert_allclose(result, [[2.0, 0.0, 0.0]])

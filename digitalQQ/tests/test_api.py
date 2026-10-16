from datetime import timedelta
import pandas as pd
import pytest
from fastapi.testclient import TestClient
from digitalQQ.backend.app import app
from digitalQQ.backend.models import Scene
from digitalQQ.backend.service import load_weather


@pytest.fixture
def client():
    with TestClient(app) as c:
        yield c


def test_api_validation_and_worker(client):
    assert client.get("/api/health").json()["status"] == "ok"
    s = Scene(magnetic=False).model_dump(mode="json")
    assert client.post("/api/scene/validate", json=s).json() == s
    result = client.post("/api/solve-pointing", json=s)
    assert result.status_code == 200
    assert result.json()["solution"]["reachable"]
    result = client.post("/api/environment?part=boundaries", json=s)
    assert result.status_code == 200
    assert result.json()["field_lines"] == []
    s["weather"]["pressure"] = 0
    assert client.post("/api/attitude", json=s).status_code == 422


def test_weather_missing_and_future(monkeypatch):
    import load

    monkeypatch.setattr(load, "load_spedas", lambda *args, **kwargs: pd.DataFrame())
    with pytest.raises(ValueError, match="没有可用"):
        load_weather("2025-01-01T00:00:00Z")
    with pytest.raises(ValueError, match="未来"):
        load_weather("2100-01-01T00:00:00Z")


def test_weather_complete_provenance(monkeypatch):
    import load

    t = pd.Timestamp("2025-01-01T00:00:00")
    monkeypatch.setattr(
        load,
        "load_spedas",
        lambda *a, **k: pd.DataFrame(
            {
                "Datetime": [t],
                "Pressure": [2.0],
                "BY_GSM": [1.0],
                "BZ_GSM": [-3.0],
                "Beta": [1.0],
                "Mgs_mach_num": [8.0],
            }
        ),
    )
    monkeypatch.setattr(
        load, "load_indices", lambda *a, **k: pd.DataFrame({"Datetime": [t], "dst": [-20.0]})
    )
    w = load_weather("2025-01-01T00:00:00Z", download=True)
    assert w["dst"] == -20 and w["source"] == "OMNI HRO + Kyoto Dst"
    assert w["observed_at"] == "2025-01-01T00:00:00+00:00"

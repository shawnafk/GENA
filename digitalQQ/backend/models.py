"""Explicit units, frames, UTC timestamps and self-contained scene inputs."""

from datetime import datetime, timezone
from typing import Literal, Annotated
import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

Vector = tuple[float, float, float]
Quaternion = tuple[float, float, float, float]


class Model(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)


class Weather(Model):
    pressure: float = Field(default=2, gt=0, le=100)
    by: float = Field(default=0, ge=-100, le=100)
    bz: float = Field(default=-3, ge=-100, le=100)
    dst: float = Field(default=-20, ge=-1000, le=200)
    beta: float = Field(default=1, ge=0, le=100)
    mach: float = Field(default=8, gt=1, le=100)
    source: str = "manual"
    observed_at: str | None = None


class OrbitSample(Model):
    time: datetime
    position: Vector
    quaternion: Quaternion | None = None

    @field_validator("time")
    @classmethod
    def utc(cls, v):
        if v.tzinfo is None:
            raise ValueError("time must include timezone, e.g. Z or +08:00")
        return v.astimezone(timezone.utc)

    @field_validator("quaternion")
    @classmethod
    def normalized(cls, v):
        if v is not None and not np.isclose(np.linalg.norm(v), 1, atol=1e-5):
            raise ValueError("quaternion must be unit xyzw, QQ -> ECI2000")
        return v


class Orbit(Model):
    source: Literal["moon", "imported"] = "moon"
    frame: Literal["ECI2000", "GSM", "GSE"] = "ECI2000"
    unit: Literal["km", "m", "RE"] = "km"
    samples: list[OrbitSample] = Field(default_factory=list, max_length=20000)

    @model_validator(mode="after")
    def coverage(self):
        if self.source == "imported":
            if len(self.samples) < 2:
                raise ValueError("imported orbit needs at least two samples")
            times = [s.time for s in self.samples]
            if any(b <= a for a, b in zip(times, times[1:])):
                raise ValueError("orbit timestamps must be strictly increasing")
            supplied = [s.quaternion is not None for s in self.samples]
            if any(supplied) and not all(supplied):
                raise ValueError("provide quaternions for every sample or none")
        return self


class Scene(Model):
    version: Literal[1] = 1
    time: datetime = datetime(2026, 10, 20, 16, tzinfo=timezone.utc)
    reference_time: datetime = datetime(2026, 10, 20, 16, tzinfo=timezone.utc)
    weather: Weather = Field(default_factory=Weather)
    orbit: Orbit = Field(default_factory=Orbit)
    roll: float = Field(default=-30.932074, ge=-180, le=180)
    pitch: float = Field(default=-1.643111, ge=-180, le=180)
    target: Vector = (-30, 0, 0)
    target_mode: Literal["point", "far", "sheet"] = "point"
    hold: Literal["earth", "inertial", "tracking"] = "earth"
    half_thickness: float = Field(default=3, gt=0, le=20)
    magnetic: bool = True
    quality: Literal["preview", "fine"] = "preview"

    @field_validator("time", "reference_time")
    @classmethod
    def utc(cls, v):
        v = OrbitSample.utc(v)
        if (
            not datetime(1950, 1, 1, tzinfo=timezone.utc)
            <= v
            <= datetime(2030, 1, 1, tzinfo=timezone.utc)
        ):
            raise ValueError(
                "场景时刻超出本地 DE432s / IGRF14 共同范围：1950-01-01 至 2030-01-01 UTC"
            )
        return v


class ExportRequest(Model):
    scene: Scene
    start: datetime
    stop: datetime
    cadence_minutes: int = Field(default=10, ge=1, le=1440)

    @field_validator("start", "stop")
    @classmethod
    def utc(cls, v):
        return OrbitSample.utc(v)

    @model_validator(mode="after")
    def interval(self):
        if not 0 < (self.stop - self.start).total_seconds() <= 7 * 86400:
            raise ValueError("export interval must be positive and at most seven days")
        if (self.stop - self.start).total_seconds() / 60 / self.cadence_minutes > 2000:
            raise ValueError("at most 2000 export samples")
        return self

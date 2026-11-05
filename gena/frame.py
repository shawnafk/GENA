"""GENA particle-event and engineering-parameter frame types."""

from __future__ import annotations

import numpy as np
from scipy.spatial.transform import Rotation

from common.frame import TimeIndexedDataFrame, TimeSeriesDataFrame
from common.time import Time
from common.interp import i_pos, i_rot

__all__ = ["ParticleDataFrame", "ParameterDataFrame"]


class _PoseAccessors:
    """Extract stored position and attitude from GENA tables without interpolation."""

    def pos(self):
        """Return row-aligned Moon-relative MJ2000 positions, (N, 3), in metres."""
        return self[["MJ2000_Star_X", "MJ2000_Star_Y", "MJ2000_Star_Z"]].to_numpy(dtype=float, copy=True)

    def rot(self):
        """Return row-aligned QQ→ECI2000 xyzw quaternions, (N, 4); invalid rows are NaN."""
        vector = self[["MJ2000_q1", "MJ2000_q2", "MJ2000_q3"]].to_numpy(dtype=float)
        squared = np.sum(vector**2, axis=1)
        # Allow float32 rounding at unit norm, as in raw attitude conversion.
        valid = np.isfinite(squared) & (squared <= 1.0 + 1e-6)
        scalar = np.sqrt(np.clip(1.0 - squared, 0.0, 1.0))
        quaternion = np.column_stack([vector, scalar])
        quaternion[~valid] = np.nan
        return quaternion


def _evaluation_timestamps(times):
    """Return one-dimensional finite Unix-second evaluation coordinates."""
    normalized = Time(times)
    values = np.atleast_1d(normalized.to_timestamp()).astype(float, copy=False)
    if not np.isfinite(values).all():
        raise ValueError("evaluation times must be finite")
    return values


class ParticleDataFrame(_PoseAccessors, TimeIndexedDataFrame):
    """GENA particle events with row-aligned position and attitude extraction."""

    def window_ave(self, key, window="5min", **kwargs):
        """Delegate time-windowed averages to :mod:`gena.post.statistic`."""
        from gena.post import statistic

        return statistic.window_ave(self, key, window=window, **kwargs)

    def window_hist0d(self, window="5min", **kwargs):
        """Delegate time-windowed event counts to :mod:`gena.post.statistic`."""
        from gena.post import statistic

        return statistic.window_hist0d(self, window=window, **kwargs)

    def window_hist1d(self, value_name, value_bins, window="5min", **kwargs):
        """Delegate time-windowed one-dimensional histograms."""
        from gena.post import statistic

        return statistic.window_hist1d(
            self, value_name, value_bins, window=window, **kwargs
        )

    def window_hist2d(self, value_names, value_bins, window="5min", **kwargs):
        """Delegate time-windowed two-dimensional histograms."""
        from gena.post import statistic

        return statistic.window_hist2d(
            self, value_names, value_bins, window=window, **kwargs
        )


class ParameterDataFrame(_PoseAccessors, TimeSeriesDataFrame):
    """GENA PDF engineering parameters with position and attitude interpolation."""

    def _loc_from_pdf(self, gap=1):
        """Return sorted valid MJ2000 position samples from this PDF frame."""
        if not isinstance(gap, (int, np.integer)) or gap < 1:
            raise ValueError("gap must be a positive integer")
        required = [
            "Timestamp",
            "MJ2000_Star_X",
            "MJ2000_Star_Y",
            "MJ2000_Star_Z",
        ]
        missing = [name for name in required if name not in self.columns]
        if missing:
            raise KeyError(f"position columns are missing: {missing}")

        values = self[required].to_numpy(dtype=float)
        valid = np.isfinite(values).all(axis=1) & (values[:, 1:] != 0).all(axis=1)
        values = values[valid]
        if len(values) == 0:
            return np.array([], dtype=float), np.empty((0, 3), dtype=float)
        values = values[np.argsort(values[:, 0], kind="stable")]
        _, unique = np.unique(values[:, 0], return_index=True)
        values = values[unique][::gap]
        return values[:, 0], values[:, 1:]

    def _rot_from_pdf(self, gap=1, norm_tol=5e-7):
        """Return sorted valid MJ2000 rotations from this PDF frame."""
        if not isinstance(gap, (int, np.integer)) or gap < 1:
            raise ValueError("gap must be a positive integer")
        if not np.isfinite(norm_tol) or norm_tol < 0:
            raise ValueError("norm_tol must be a finite non-negative number")
        required = ["Timestamp", "MJ2000_q1", "MJ2000_q2", "MJ2000_q3"]
        missing = [name for name in required if name not in self.columns]
        if missing:
            raise KeyError(f"rotation columns are missing: {missing}")

        values = self[required].to_numpy(dtype=float)
        vector = values[:, 1:]
        squared = np.sum(vector**2, axis=1)
        scalar = np.sqrt(1.0 - np.clip(squared, 0.0, 1.0))
        quaternion = np.column_stack([vector, scalar])
        norms = np.linalg.norm(quaternion, axis=1)
        valid = np.isfinite(values).all(axis=1)
        present = (vector != 0).any(axis=1)
        if "PosePresent" in self.columns:
            # Legacy PDFs omit the flag and encode absent attitude as three zeros.
            flag = self["PosePresent"].astype("boolean")
            present = np.where(flag.isna(), present, flag.to_numpy(dtype=bool, na_value=False))
        valid &= present
        valid &= np.isfinite(norms)
        valid &= (norms >= 1.0 - norm_tol) & (norms <= 1.0 + norm_tol)
        timestamps = values[valid, 0]
        quaternion = quaternion[valid]
        if len(timestamps) == 0:
            raise ValueError("no valid PDF rotation samples")

        quaternion = quaternion / np.linalg.norm(quaternion, axis=1, keepdims=True)
        order = np.argsort(timestamps, kind="stable")
        timestamps = timestamps[order]
        quaternion = quaternion[order]
        _, unique = np.unique(timestamps, return_index=True)
        timestamps = timestamps[unique][::gap]
        quaternion = quaternion[unique][::gap]
        return timestamps.tolist(), Rotation.from_quat(quaternion)

    def pos_at(self, times, *, bounds_error=True):
        """Linearly interpolate valid PDF positions at arbitrary UTC times."""
        evaluation = _evaluation_timestamps(times)
        if len(evaluation) == 0:
            return np.empty((0, 3), dtype=float)
        source, positions = self._loc_from_pdf()
        if len(source) == 0:
            raise ValueError("no valid PDF position samples")
        covered = (evaluation >= source[0]) & (evaluation <= source[-1])
        if bounds_error and not covered.all():
            if len(source) == 1:
                raise ValueError("position interpolation requires at least two source samples")
            raise ValueError("position evaluation times are outside PDF coverage")
        result = np.full((len(evaluation), 3), np.nan)
        if len(source) == 1:
            result[covered] = positions
        else:
            result[covered] = i_pos(source, positions, evaluation[covered], extrapolate=False)
        return result

    def rot_at(self, times):
        """SLERP valid PDF attitudes at arbitrary UTC times."""
        evaluation = _evaluation_timestamps(times)
        if len(evaluation) == 0:
            raise ValueError("rotation evaluation times must contain at least one value")
        source, rotations = self._rot_from_pdf()
        source = np.asarray(source, dtype=float)
        if np.any(evaluation < source[0]) or np.any(evaluation > source[-1]):
            raise ValueError("rotation evaluation times are outside PDF coverage")
        order = np.argsort(evaluation, kind="stable")
        interpolated = i_rot(source, rotations, evaluation[order])
        inverse = np.empty_like(order)
        inverse[order] = np.arange(len(order))
        return interpolated[inverse]

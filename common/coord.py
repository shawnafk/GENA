"""Cartesian coordinates, vector geometry and rotation helpers.

Coordinates stores only NumPy arrays and metadata. Standard frame conversions
use lazy SpacePy, Astropy or PySPEDAS adapters; QQ/A/B and fixed-angle AGSM
are shared project transforms. Conversion never changes the reference origin.

QQ looks along -Z. Probe-to-QQ mounting rotations are -5 degrees for A and
+5 degrees for B around X. AGSM uses a fixed 4.2-degree GSM XY rotation.
This module also owns FOV, angle, quaternion application and SLERP operations.
"""

from __future__ import annotations

import numpy as np

from common.time import Time
from common._coord_backends import BACKENDS, transform as _backend_transform, validate_backend

#very related constants, so directly defined here rather than in common.constants
AGSM_ANGLE_DEG = 4.2
_AGSM_ANGLE_RAD = np.radians(AGSM_ANGLE_DEG)
_R_GSM_TO_AGSM = np.array([
    [np.cos(_AGSM_ANGLE_RAD), -np.sin(_AGSM_ANGLE_RAD), 0.0],
    [np.sin(_AGSM_ANGLE_RAD), np.cos(_AGSM_ANGLE_RAD), 0.0],
    [0.0, 0.0, 1.0],
])
_R_AGSM_TO_GSM = _R_GSM_TO_AGSM.T


def _vectors(values):
    """Validate a single vector or row-wise vector batch without changing rank."""
    array = np.asarray(values, dtype=float)
    if array.shape == (3,) or (array.ndim == 2 and array.shape[1] == 3):
        return array
    raise ValueError(f"invalid input vector shape: {array.shape}, expected (3,) or (N, 3)")


def vector_to_fov(vectors, axes=None, *, direction='propagation'):
    """Return px/py angles in degrees for basis rows and a selected boresight.

    Input is (3,) or (N, 3); basis is (3, 3) or sample-wise (N, 3, 3).
    With no basis, propagation uses +Z and look uses -Z, with +X/+Y lateral
    axes. Explicit basis rows define the convention and cannot be combined
    with direction='look'. Zero/non-finite vectors are rejected.
    """
    vectors = _vectors(vectors)
    if direction not in {'propagation', 'look'}:
        raise ValueError("direction must be 'propagation' or 'look'")
    if not np.isfinite(vectors).all() or np.any(np.linalg.norm(vectors, axis=-1) == 0):
        raise ValueError("FOV vectors must be finite and nonzero")
    if axes is None:
        basis = np.diag([1., 1., -1. if direction == 'look' else 1.])
    else:
        if direction != 'propagation':
            raise ValueError("explicit axes already define the FOV convention")
        basis = np.asarray(axes, dtype=float)
    if basis.shape == (3, 3):
        projected = vectors @ basis.T
    elif vectors.ndim == 2 and basis.shape == (len(vectors), 3, 3):
        projected = np.einsum('nij,nj->ni', basis, vectors)
    else:
        raise ValueError("axes must have shape (3, 3) or (N, 3, 3)")
    if not np.isfinite(basis).all() or not np.allclose(basis @ np.swapaxes(basis, -1, -2), np.eye(3)):
        raise ValueError("FOV basis rows must be orthonormal")
    return (np.degrees(np.arctan2(projected[..., 0], projected[..., 2])),
            np.degrees(np.arctan2(projected[..., 1], projected[..., 2])))


def fov_to_vector(px, py, *, direction='propagation'):
    """Return unit vectors for degree angles strictly inside (-90, 90).

    Scalar angles return (3,); broadcast one-dimensional inputs return (N, 3).
    Zero angles point along +Z for propagation or -Z for look. The inverse
    is vector_to_fov with the same direction, within the forward hemisphere.
    """
    if direction not in {'propagation', 'look'}:
        raise ValueError("direction must be 'propagation' or 'look'")
    px, py = np.broadcast_arrays(np.asarray(px, dtype=float), np.asarray(py, dtype=float))
    if px.ndim > 1 or not np.isfinite(px).all() or not np.isfinite(py).all():
        raise ValueError("FOV angles must be finite scalars or one-dimensional arrays")
    if np.any(np.abs(px) >= 90) or np.any(np.abs(py) >= 90):
        raise ValueError("FOV angles must lie strictly inside (-90, 90) degrees")
    vectors = np.stack([np.tan(np.radians(px)), np.tan(np.radians(py)),
                        np.full(px.shape, -1. if direction == 'look' else 1.)], axis=-1)
    return vectors / np.linalg.norm(vectors, axis=-1, keepdims=True)


def angle_between(a, b, *, degrees=True):
    """Return the unsigned angle between broadcastable nonzero vectors."""
    a, b = np.broadcast_arrays(_vectors(a), _vectors(b))
    denominator = np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1)
    if not np.isfinite(a).all() or not np.isfinite(b).all() or np.any(denominator == 0):
        raise ValueError("angle vectors must be finite and nonzero")
    cosine = np.sum(a * b, axis=-1) / denominator
    result = np.arccos(np.clip(cosine, -1., 1.))
    return np.degrees(result) if degrees else result


def angle_to_plane(vectors, normal, *, degrees=True, signed=True):
    """Return angle to a plane; positive sign follows its oriented normal.

    Signed results lie in [-90, 90] degrees; signed=False returns magnitudes.
    Plane offsets do not affect a vector-to-plane angle.
    """
    angle = np.pi / 2 - angle_between(vectors, normal, degrees=False)
    if not signed:
        angle = np.abs(angle)
    return np.degrees(angle) if degrees else angle


def quaternion_array(rotation, *, count=None):
    """Return copied unit xyzw quaternion rows from arrays or SciPy Rotation.

    Accept (4,), (1, 4), or (N, 4). An optional count checks whether rotations
    are single/broadcastable or sample-aligned. No frame meaning is assigned.
    """
    from scipy.spatial.transform import Rotation

    raw = rotation.as_quat() if isinstance(rotation, Rotation) else np.asarray(rotation, dtype=float)
    if raw.shape == (4,):
        raw = raw[None, :]
    if raw.ndim != 2 or raw.shape[1] != 4:
        raise ValueError('quaternions must have shape (4,), (1, 4), or (N, 4)')
    if count is not None and len(raw) not in {1, count}:
        raise ValueError('quaternion count must be one or match the sample count')
    norm = np.linalg.norm(raw, axis=1, keepdims=True)
    if not np.isfinite(raw).all() or np.any(norm == 0):
        raise ValueError('quaternions must be finite and nonzero')
    return raw / norm


def rotate_vectors(vectors, rotation, *, inverse=False):
    """Apply a Cartesian rotation, preserving (3,) or (N, 3) vector shape.

    rotation is a fixed (3, 3) matrix, xyzw quaternion(s), or SciPy Rotation.
    Matrix inputs use the column-vector convention R @ v; inverse uses R.T.
    Quaternion inputs may be single or sample-aligned. No translation occurs.
    """
    vectors = _vectors(vectors)
    if isinstance(rotation, (np.ndarray, list, tuple)) and np.shape(rotation) == (3, 3):
        matrix = np.asarray(rotation, dtype=float)
        return vectors @ (matrix if inverse else matrix.T)
    from scipy.spatial.transform import Rotation

    quaternion = quaternion_array(rotation, count=len(np.atleast_2d(vectors)))
    if not vectors.size:
        return vectors.copy()
    if vectors.ndim == 1:
        quaternion = quaternion[0]
    operator = Rotation.from_quat(quaternion)
    return (operator.inv() if inverse else operator).apply(vectors)


def interpolate_quaternions(source, rotation, target):
    """SLERP xyzw rotations on strictly increasing numeric source coordinates.

    Axes are numeric interpolation coordinates in a common unit, not parsed
    dates. Targets may be unordered or repeated but must be inside coverage.
    A single rotation is constant; source sample cleanup belongs to the caller.
    Returns (len(target), 4), retaining no SciPy objects.
    """
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    if source.ndim != 1 or target.ndim != 1 or not len(source):
        raise ValueError('source and target must be 1-D, with nonempty source')
    if (not np.isfinite(source).all() or not np.isfinite(target).all()
            or np.any(np.diff(source) <= 0)):
        raise ValueError('source must be finite and strictly increasing; targets must be finite')
    if np.any(target < source[0]) or np.any(target > source[-1]):
        raise ValueError('rotation target coordinates are outside source coverage')
    quaternion = quaternion_array(rotation, count=len(source))
    if len(quaternion) == 1:
        return np.repeat(quaternion, len(target), axis=0)
    if not len(target):
        return np.empty((0, 4))
    from scipy.spatial.transform import Rotation, Slerp

    return Slerp(source, Rotation.from_quat(quaternion))(target).as_quat()


def _times(times, count=None):
    """Normalize UTC times, broadcasting only a scalar time when count is given."""
    time = Time(times)
    if time.isscalar:
        dates = np.full(1 if count is None else count, time.datetime.to_datetime64(),
                        dtype='datetime64[ns]')
    else:
        dates = time.datetime.to_numpy(dtype='datetime64[ns]', copy=True)
    if count is not None and len(dates) != count:
        raise ValueError("times and values must have the same sample count")
    return dates


def _probe_to_qq(values, probe, *, inverse=False):
    """Apply the fixed ±5-degree probe mounting rotation in either direction."""
    angle = np.radians(-5.0 if probe == 'A' else 5.0)
    if inverse:
        angle = -angle
    cosine, sine = np.cos(angle), np.sin(angle)
    matrix = np.array([[1., 0., 0.], [0., cosine, -sine], [0., sine, cosine]])
    return rotate_vectors(values, matrix)

def transform_vectors(values, *, source, target, times=None, backend='spacepy',
                      attitude=None, unit=None, origin='earth'):
    """Transform Cartesian arrays through one public conversion entry.

    Parameters
    ----------
    values : array-like, shape (3,) or (N, 3)
        Components. Result preserves single-vector or row-wise batch shape.
    source, target : str
        Canonical uppercase frame names (e.g. ECI2000, GSM), used as supplied.
        QQ/A/B and AGSM are project-defined; standard
        frame support is determined by the explicitly selected backend.
    times : scalar or one-dimensional time-like, optional
        UTC times required only when a standard backend transform is needed.
        Scalars broadcast; sequences must match the input sample count.
    backend : {'spacepy', 'astropy', 'pyspedas'}
        Explicit backend; unsupported conversions never fall back.
    attitude : array-like or scipy Rotation, optional
        QQ→ECI2000 xyzw quaternion(s), required when crossing this boundary.
    unit : str or None
        Preserved component unit; Astropy position transforms need length units.
    origin : str
        Preserved reference origin, default earth; no implicit translation.

    Returns
    -------
    numpy.ndarray
        New component array in the target frame, with the original input shape.
        Local A/B/QQ and GSM/AGSM rotations do not need times or backend packages.
    """
    backend = validate_backend(backend)
    original = _vectors(values)
    result = np.atleast_2d(original).copy()
    dates = None if times is None else _times(times, len(result))
    quat = None if attitude is None else quaternion_array(attitude, count=len(result))
    local = {'QQ', 'A', 'B'}
    if source != target:
        if source in local and target in local:
            if source != 'QQ':
                result = _probe_to_qq(result, source)
            if target != 'QQ':
                result = _probe_to_qq(result, target, inverse=True)
        else:
            if source in local:
                if source != 'QQ':
                    result = _probe_to_qq(result, source)
                if quat is None:
                    raise ValueError('QQ/A/B ↔ standard frames requires attitude (QQ → ECI2000)')
                # Validated CE7 convention: QQ→ECI2000 is R @ v.
                result = rotate_vectors(result, quat)
                source = 'ECI2000'
            elif source == 'AGSM':
                result = rotate_vectors(result, _R_AGSM_TO_GSM)
                source = 'GSM'
            relay = 'ECI2000' if target in local else 'GSM' if target == 'AGSM' else target
            if source != relay:
                if dates is None:
                    raise ValueError('times are required for standard-frame conversion')
                result = _backend_transform(result, dates, source=source, target=relay,
                                            unit=unit, origin=str(origin).lower(), backend=backend)
            if target in local:
                if quat is None:
                    raise ValueError('QQ/A/B ↔ standard frames requires attitude (QQ → ECI2000)')
                result = rotate_vectors(result, quat, inverse=True)
                if target != 'QQ':
                    result = _probe_to_qq(result, target, inverse=True)
            elif target == 'AGSM':
                result = rotate_vectors(result, _R_GSM_TO_AGSM)
    return result[0] if original.ndim == 1 else result


class Coordinates:
    """Cartesian samples stored independently of any coordinate backend.

    Parameters
    ----------
    values : array-like, shape (3,) or (N, 3)
        Cartesian components, stored as a copied float array of shape (N, 3).
    times : scalar or one-dimensional time-like
        UTC times normalized by common.time.Time to NumPy datetime64[ns].
        A scalar is broadcast; sequences must match the sample count.
    frame : str
        Source frame. QQ/A/B/AGSM are project frames; standard support depends
        on the backend chosen when calling to(). Use canonical uppercase names
        such as ECI2000 and GSM; no spelling or alias normalization is applied.
    unit : str or None
        Numeric component unit, preserved by conversion. Astropy native
        position transforms require an explicit length unit (including Re).
    origin : str
        Reference origin, default earth. Conversion rotates/transforms within
        this origin; it never silently adds spacecraft or lunar offsets.
    attitude : array-like or scipy Rotation, optional
        QQ→ECI2000 xyzw quaternion(s), stored only as a NumPy array. Needed
        when crossing between QQ/A/B and standard frames.

    No backend objects or converted-frame cache are stored. to() and
    resample_at() return new objects. Inputs may contain missing vector values.
    """

    def __init__(self, values, times, frame, *, unit=None, origin='earth', attitude=None):
        self.values = np.atleast_2d(_vectors(values)).copy()
        self.times = _times(times, len(self.values))
        self.frame = frame
        self.unit = unit
        self.origin = str(origin).lower()
        self.attitude = None if attitude is None else quaternion_array(attitude, count=len(self.values))

    def __len__(self):
        return len(self.values)

    def __repr__(self):
        return f"Coordinates(N={len(self)}, frame={self.frame!r}, unit={self.unit!r}, origin={self.origin!r})"

    def to(self, frame, *, backend='spacepy', attitude=None):
        """Return converted coordinates using an explicit backend, without fallback.

        SpacePy and PySPEDAS support space-physics frames. Astropy supports
        native ITRS/GCRS/CIRS/TEME and ICRS/Galactic positions, but does not
        pretend GCRS is ECI2000 or provide a native GSM/GSE alias. Local probe
        and AGSM rotations do not require an installed backend package.
        """
        target = frame
        quat = self.attitude if attitude is None else quaternion_array(attitude, count=len(self))
        values = transform_vectors(self.values, times=self.times, source=self.frame,
                                   target=target, backend=backend, attitude=quat,
                                   unit=self.unit, origin=self.origin)
        return Coordinates(values, self.times, target, unit=self.unit,
                           origin=self.origin, attitude=quat)

    def resample_at(self, times):
        """Linearly sample native-frame components and SLERP attitude, without extrapolation.

        Sort source times and keep the first duplicate. Preserve target order
        and duplicates. Targets outside coverage raise ValueError. Interpolation
        is of components in this frame, not of already converted frame caches.
        """
        dates = _times(times)
        if not len(dates):
            quat = None if self.attitude is None else self.attitude[:1]
            return Coordinates(np.empty((0, 3)), dates, self.frame, unit=self.unit,
                               origin=self.origin, attitude=quat)
        if not len(self):
            raise ValueError("cannot resample an empty coordinate series")
        order = np.argsort(self.times, kind='stable')
        source, first = np.unique(self.times[order], return_index=True)
        indices = order[first]
        if np.any(dates < source[0]) or np.any(dates > source[-1]):
            raise ValueError("target times are outside source coverage")
        x = (source - source[0]) / np.timedelta64(1, 's')
        target = (dates - source[0]) / np.timedelta64(1, 's')
        values = np.column_stack([np.interp(target, x, self.values[indices, i]) for i in range(3)])
        quat = self.attitude
        if quat is not None:
            if len(quat) > 1:
                quat = quat[indices]
            quat = interpolate_quaternions(x, quat, target)
        return Coordinates(values, dates, self.frame, unit=self.unit,
                           origin=self.origin, attitude=quat)

    def _combine(self, other, operation):
        """Combine aligned components only when frame, origin and units agree."""
        if not isinstance(other, Coordinates):
            return NotImplemented
        if (self.frame, self.unit, self.origin) != (other.frame, other.unit, other.origin):
            raise ValueError("coordinate arithmetic requires matching frame, unit and origin")
        if not np.array_equal(self.times, other.times):
            raise ValueError("coordinate arithmetic requires identical time axes")
        if self.frame in {'QQ', 'A', 'B'}:
            left = None if self.attitude is None else np.broadcast_to(self.attitude, (len(self), 4))
            right = None if other.attitude is None else np.broadcast_to(other.attitude, (len(other), 4))
            if (left is None) != (right is None) or (left is not None and not np.allclose(
                    np.abs(np.sum(left * right, axis=1)), 1.0)):
                raise ValueError("local-frame arithmetic requires matching attitudes")
        return Coordinates(operation(self.values, other.values), self.times, self.frame,
                           unit=self.unit, origin=self.origin, attitude=self.attitude)

    def __add__(self, other):
        return self._combine(other, np.add)

    def __sub__(self, other):
        return self._combine(other, np.subtract)


class coord(Coordinates):
    """Compatibility entry for coord(times, values, frame, rot).

    New code should use Coordinates(values, times, frame). convert() fills the
    legacy array dictionary; standard conversions default to SpacePy. No
    backend object is retained. Legacy absent attitude means identity.
    """

    def __init__(self, dt, coord_value, coord_type, rot=None, *, backend='spacepy', unit=None, origin='earth'):
        self.backend = validate_backend(backend)
        if isinstance(coord_value, Coordinates):
            coord_value = coord_value.to(coord_type, backend=backend).values
        if rot is None or (np.isscalar(rot) and rot == 0):
            rot = np.array([0., 0., 0., 1.])
        super().__init__(coord_value, dt, coord_type, unit=unit, origin=origin, attitude=rot)
        self.c = {self.frame: self.values.copy()}

    @property
    def datetime(self):
        """Legacy time view; canonical storage is the NumPy times array."""
        return Time(self.times).datetime

    def convert(self, out_type, *, backend=None):
        """Fill a legacy array result, recomputing from the canonical samples."""
        target = out_type
        result = self.to(target, backend=self.backend if backend is None else backend)
        self.c[target] = result.values

    def sync(self, new_datetime):
        """Replace canonical samples by resampled data and discard old conversions."""
        result = self.resample_at(new_datetime)
        self.values, self.times, self.attitude = result.values, result.times, result.attitude
        self.c = {self.frame: self.values.copy()}

    def __getitem__(self, key):
        if key not in self.c:
            self.convert(key)
        return self.c[key]

    def _combine(self, other, operation):
        result = super()._combine(other, operation)
        if result is NotImplemented:
            return result
        return coord(result.times, result.values, result.frame, result.attitude,
                     backend=self.backend, unit=result.unit, origin=result.origin)


__all__ = [
    'Coordinates', 'coord', 'transform_vectors',
    'angle_between', 'angle_to_plane', 'vector_to_fov', 'fov_to_vector',
    'quaternion_array', 'rotate_vectors', 'interpolate_quaternions',
]

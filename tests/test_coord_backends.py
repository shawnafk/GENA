"""Coordinate data boundaries, actual adapters and legacy interoperability."""

import subprocess
import sys

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from common.coord import Coordinates, coord, transform_vectors

TIMES = ['2024-01-01T00:00:00', '2024-01-01T00:00:10']
VALUES = np.array([[2., 1., -3.], [4., -2., 1.]])


def test_import_and_local_operations_do_not_import_backends():
    code = '''
import sys
from common.time import Time
baseline = set(sys.modules)
from common.coord import Coordinates
x = Coordinates([1, 2, 3], '2024-01-01', 'GSM')
x.to('AGSM', backend='astropy')
assert not any(k in set(sys.modules) - baseline for k in ('spacepy', 'astropy.coordinates', 'pyspedas', 'libs.utils'))
assert not hasattr(x, 'temp')
'''
    subprocess.run([sys.executable, '-c', code], check=True)


def test_array_storage_normalizes_time_and_copies_input():
    raw = VALUES.copy()
    x = Coordinates(raw, '2024-01-01T08:00:00+08:00', 'GSM')
    raw[:] = 0
    assert x.values == pytest.approx(VALUES)
    assert x.times.dtype == np.dtype('datetime64[ns]')
    assert np.all(x.times == np.datetime64('2024-01-01T00:00:00'))
    with pytest.raises(ValueError, match='sample count'):
        Coordinates(VALUES, [TIMES[0]], 'GSM')
    with pytest.raises(ValueError):
        Coordinates([1, 2, 3], [None], 'GSM')


@pytest.mark.parametrize('backend', ['spacepy', 'pyspedas'])
def test_actual_space_physics_backend_round_trip(backend):
    pytest.importorskip(backend)
    x = Coordinates(VALUES, TIMES, 'ECI2000', unit='Re')
    gsm = x.to('GSM', backend=backend)
    restored = gsm.to('ECI2000', backend=backend)
    assert restored.values == pytest.approx(VALUES, abs=1e-7)
    assert gsm.unit == 'Re'
    assert not hasattr(gsm, 'backend')
    assert x.values == pytest.approx(VALUES)


def test_astropy_native_positions_round_trip():
    pytest.importorskip('astropy')
    x = Coordinates(VALUES * 1000, TIMES, 'ITRS', unit='km')
    y = x.to('GCRS', backend='astropy')
    assert y.to('ITRS', backend='astropy').values == pytest.approx(x.values, abs=1e-7)
    with pytest.raises(NotImplementedError, match='native'):
        x.to('GSM', backend='astropy')
    with pytest.raises(ValueError, match='length unit'):
        Coordinates(VALUES, TIMES, 'ITRS').to('GCRS', backend='astropy')


def test_astropy_re_unit_and_origin_checks():
    pytest.importorskip('astropy')
    x = Coordinates(VALUES, TIMES, 'GEO', unit='Re')
    assert x.to('GCRS', backend='astropy').to('GEO', backend='astropy').values == pytest.approx(VALUES)
    with pytest.raises(NotImplementedError, match='origins'):
        x.to('ICRS', backend='astropy')
    with pytest.raises(NotImplementedError, match='origin'):
        Coordinates(VALUES, TIMES, 'GSM', origin='moon').to('GSE')


@pytest.mark.parametrize('probe', ['A', 'B', 'QQ'])
@pytest.mark.parametrize('batch', [False, True])
def test_attitude_single_and_batch_round_trip(probe, batch):
    attitude = Rotation.from_euler('z', [20., 40.] if batch else 30., degrees=True)
    x = Coordinates(VALUES, TIMES, probe, attitude=attitude)
    expected = transform_vectors(VALUES, source=probe, target='QQ') if probe != 'QQ' else VALUES
    inertial = x.to('ECI2000')
    assert inertial.values == pytest.approx(attitude.apply(expected))
    assert inertial.to(probe).values == pytest.approx(VALUES)


def test_attitude_required_only_for_crossing_inertial_boundary():
    x = Coordinates(VALUES, TIMES, 'A')
    assert x.to('B').to('A').values == pytest.approx(VALUES)
    with pytest.raises(ValueError, match='requires attitude'):
        x.to('ECI2000')
    with pytest.raises(ValueError, match='backend'):
        x.to('A', backend='unknown')


def test_resampling_slerp_and_legacy_sync_use_new_canonical_data():
    rotation = Rotation.from_euler('z', [0., 90.], degrees=True)
    x = Coordinates(VALUES, TIMES, 'QQ', attitude=rotation)
    sampled = x.resample_at(['2024-01-01T00:00:05'])
    expected = Rotation.from_euler('z', 45., degrees=True).apply(VALUES.mean(axis=0))
    assert sampled.to('ECI2000').values[0] == pytest.approx(expected)
    assert len(x) == 2
    legacy = coord(TIMES, VALUES, 'QQ', rotation)
    legacy.convert('ECI2000')
    legacy.sync(['2024-01-01T00:00:05'])
    assert set(legacy.c) == {'QQ'}
    assert legacy['ECI2000'][0] == pytest.approx(expected)
    with pytest.raises(ValueError, match='coverage'):
        x.resample_at(['2023-01-01'])


def test_resampling_unsorted_duplicates_empty_and_singleton():
    x = Coordinates([[4, 0, 0], [0, 0, 0], [99, 0, 0]],
                    [TIMES[1], TIMES[0], TIMES[0]], 'GSM')
    result = x.resample_at([TIMES[1], TIMES[0], TIMES[1]])
    assert result.values[:, 0] == pytest.approx([4, 0, 4])
    assert x.resample_at([]).values.shape == (0, 3)
    assert Coordinates([1, 2, 3], TIMES[0], 'GSM').resample_at([TIMES[0]]).values.shape == (1, 3)


def test_arithmetic_rejects_unaligned_metadata_and_times():
    x = Coordinates(VALUES, TIMES, 'GSM', unit='Re')
    assert (x - x).values == pytest.approx(np.zeros_like(VALUES))
    for other in [Coordinates(VALUES, TIMES, 'GSE', unit='Re'),
                  Coordinates(VALUES, TIMES, 'GSM', unit='km'),
                  Coordinates(VALUES, TIMES[::-1], 'GSM', unit='Re')]:
        with pytest.raises(ValueError):
            x + other


def test_missing_backend_does_not_prevent_storage_or_local_transforms():
    code = """
import sys
from common.coord import Coordinates
class BlockBackends:
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'spacepy', 'astropy', 'pyspedas'}:
            raise ModuleNotFoundError('blocked optional dependency', name=fullname)
sys.meta_path.insert(0, BlockBackends())
x = Coordinates([1, 2, 3], '2024-01-01', 'GSM')
x.to('AGSM')
for backend in ('spacepy', 'pyspedas', 'astropy'):
    source = Coordinates([1, 2, 3], '2024-01-01',
                         'ITRS' if backend == 'astropy' else 'GSM', unit='km')
    try:
        source.to('GCRS' if backend == 'astropy' else 'GSE', backend=backend)
    except ImportError as error:
        assert backend in str(error)
    else:
        raise AssertionError('missing dependency was not reported')
"""
    subprocess.run([sys.executable, '-c', code], check=True)


@pytest.mark.parametrize('source,target', [('GSM', 'GSE'), ('ECI2000', 'GSM')])
@pytest.mark.parametrize('backend', ['spacepy', 'pyspedas'])
def test_array_entry_and_class_share_actual_backend(source, target, backend):
    pytest.importorskip(backend)
    values = transform_vectors(VALUES, times=TIMES, source=source, target=target, backend=backend)
    assert values == pytest.approx(Coordinates(VALUES, TIMES, source).to(target, backend=backend).values)
    with pytest.raises(ValueError, match='times are required'):
        transform_vectors(VALUES, source=source, target=target, backend=backend)


@pytest.mark.parametrize('source,target', [('A', 'QQ'), ('B', 'A'), ('GSM', 'AGSM')])
def test_array_entry_preserves_rank_and_needs_no_time_for_local_frames(source, target):
    single = transform_vectors(VALUES[0], source=source, target=target)
    batch = transform_vectors(VALUES[:1], source=source, target=target)
    assert single.shape == (3,)
    assert batch.shape == (1, 3)
    assert batch[0] == pytest.approx(single)
    assert transform_vectors(single, source=target, target=source) == pytest.approx(VALUES[0])


def test_array_entry_astropy_positions_and_attitude_boundary():
    pytest.importorskip('astropy')
    values = transform_vectors(VALUES, times=TIMES, source='ITRS', target='GCRS',
                               backend='astropy', unit='km')
    assert values == pytest.approx(Coordinates(VALUES, TIMES, 'ITRS', unit='km').to('GCRS', backend='astropy').values)
    quaternion = Rotation.from_euler('y', 30., degrees=True)
    inertial = transform_vectors(VALUES[0], source='QQ', target='ECI2000', attitude=quaternion)
    assert inertial == pytest.approx(quaternion.apply(VALUES[0]))
    with pytest.raises(ValueError, match='requires attitude'):
        transform_vectors(VALUES[0], source='QQ', target='ECI2000')

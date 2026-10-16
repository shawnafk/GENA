import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import pytest

from plots.basic import themis as vis
from libs import plasma


def _products():
    esa_time = pd.date_range('2024-01-01', periods=3, freq='min')
    sst_time = pd.date_range('2024-01-01 00:00:30', periods=2, freq='2min')
    esa = {'tha_peif_en_eflux': {
        'x': esa_time,
        'v': np.array([[10., 100.], [12., 120.], [14., 140.]]),
        'y': np.array([[1., 10.], [2., 20.], [3., 30.]]),
    }}
    sst = {'tha_psir_en_eflux': {
        'x': sst_time,
        'v': np.array([30000., 100000., 200000.]),
        'y': np.array([[100., 200., 300.], [150., 250., 350.]]),
    }}
    return esa, sst


def test_esa_sst_spectra_keeps_independent_grids_and_colorbars():
    esa, sst = _products()
    figure, axes, meshes = vis.esa_sst_spectra(esa, sst, prefix='tha', sst_key='psir_en_eflux')
    assert axes[0].get_shared_x_axes().joined(*axes)
    assert meshes[0].get_array().size == 6
    assert meshes[1].get_array().size == 6
    assert all(isinstance(mesh.norm, LogNorm) for mesh in meshes)
    np.testing.assert_allclose(
        np.sort(np.asarray(meshes[0].get_array()).ravel()),
        [100., 150., 200., 250., 300., 350.],
    )
    np.testing.assert_allclose(
        np.sort(np.asarray(meshes[1].get_array()).ravel()),
        [1., 2., 3., 10., 20., 30.],
    )
    # Two panels plus one top colorbar for each panel.
    assert len(figure.axes) == 4
    plt.close(figure)


def test_shared_color_scale_and_time_averages_include_sorted_combination():
    esa, sst = _products()
    figure, _, meshes = vis.esa_sst_spectra(
        esa, sst, prefix='tha', sst_key='psir_en_eflux', shared_color_scale=True
    )
    assert meshes[0].get_clim() == meshes[1].get_clim()
    assert meshes[0].get_clim() == (1.0, 350.0)
    energy, flux = vis.time_averaged_spectra(
        esa, sst, prefix='tha', sst_key='psir_en_eflux'
    )
    np.testing.assert_allclose(
        energy,
        [12., 120., 30000., 100000., 200000.],
    )
    np.testing.assert_allclose(
        flux,
        [2., 20., 125., 225., 325.],
    )
    plt.close(figure)


def test_time_average_removes_nonfinite_energy_flux_pairs():
    esa, sst = _products()
    esa['tha_peif_en_eflux']['v'] = np.array([
        [10., 100., 1000.],
        [12., 120., np.nan],
        [14., 140., np.nan],
    ])
    esa['tha_peif_en_eflux']['y'] = np.array([
        [np.nan, 10., np.nan],
        [np.nan, 20., np.nan],
        [np.nan, 30., np.nan],
    ])

    energy, flux = vis.time_averaged_spectra(
        esa, sst, prefix='tha', sst_key='psir_en_eflux'
    )

    np.testing.assert_allclose(
        energy,
        [120., 30000., 100000., 200000.],
    )
    np.testing.assert_allclose(
        flux,
        [20., 125., 225., 325.],
    )
    assert np.isfinite(energy).all()
    assert np.isfinite(flux).all()
    assert np.all(np.diff(energy) >= 0)


def test_particle_flux_conversion_happens_before_time_average():
    esa, sst = _products()

    energy, particle_flux = vis.time_averaged_spectra(
        esa, sst, prefix='tha', sst_key='psir_en_eflux', en=False
    )

    expected_esa_bin_1 = np.mean([1.0 / 10.0, 2.0 / 12.0, 3.0 / 14.0])
    expected_esa_bin_2 = np.mean([10.0 / 100.0, 20.0 / 120.0, 30.0 / 140.0])
    np.testing.assert_allclose(
        energy,
        [12.0, 120.0, 30000.0, 100000.0, 200000.0],
    )
    np.testing.assert_allclose(
        particle_flux,
        [
            expected_esa_bin_1,
            expected_esa_bin_2,
            125.0 / 30000.0,
            225.0 / 100000.0,
            325.0 / 200000.0,
        ],
    )


def test_time_averaged_spectra_normalizes_in_native_ev_units():
    esa, sst = _products()

    energy, flux = vis.time_averaged_spectra(
        esa,
        sst,
        prefix='tha',
        sst_key='psir_en_eflux',
        normalized=True,
        E_norm=(20000.0, 40000.0),
    )

    np.testing.assert_allclose(energy, [12.0, 120.0, 30000.0, 100000.0, 200000.0])
    np.testing.assert_allclose(flux, [2.0 / 125.0, 20.0 / 125.0, 1.0, 225.0 / 125.0, 325.0 / 125.0])


def test_particle_flux_spectrogram_and_colorbar_label():
    esa, sst = _products()

    figure, _, meshes = vis.esa_sst_spectra(
        esa, sst, prefix='tha', sst_key='psir_en_eflux', en=False
    )

    assert meshes[0].get_array().size == 6
    assert meshes[1].get_array().size == 6
    colorbar_labels = {
        label
        for axis in figure.axes[2:]
        for label in (axis.get_xlabel(), axis.get_ylabel())
        if label
    }
    assert colorbar_labels == {r"$j$"}
    plt.close(figure)


def test_spectrum_en_must_be_boolean():
    esa, _ = _products()

    with pytest.raises(TypeError, match='en must be bool'):
        vis.spectrum(esa, 'peif_en_eflux', prefix='tha', en='no')


def test_time_averaged_particle_moments_preserve_density_weighted_flux():
    times = pd.date_range('2025-01-01', periods=3, freq='min')
    esa = {
        'tha_peir_density': {
            'x': times,
            'y': np.array([1.0, 2.0, 1.0]),
        },
        'tha_peir_velocity_gsm': {
            'x': times,
            'y': np.array([
                [100.0, 0.0, 0.0],
                [300.0, 0.0, 0.0],
                [100.0, 0.0, 0.0],
            ]),
        },
        'tha_peir_avgtemp': {
            'x': times,
            'y': np.full(3, 1000.0),
        },
    }

    moments = vis.time_averaged_particle_moments(esa, prefix='tha')

    assert moments['density_cm3'] == pytest.approx(4.0 / 3.0)
    np.testing.assert_allclose(
        moments['velocity_gsm_kms'], [200.0, 0.0, 0.0]
    )
    assert moments['speed_kms'] == pytest.approx(200.0)
    assert moments['temperature_eV'] > 1000.0
    assert moments['sample_count'] == 3


def test_combined_esa_sst_moments_add_partial_moments():
    times = pd.date_range('2025-01-01', periods=3, freq='min')
    esa = {
        'tha_peir_density': {'x': times, 'y': np.full(3, 1.0)},
        'tha_peir_velocity_gsm': {
            'x': times,
            'y': np.tile([100.0, 0.0, 0.0], (3, 1)),
        },
        'tha_peir_avgtemp': {'x': times, 'y': np.full(3, 1000.0)},
    }
    sst = {
        'tha_psif_density': {'x': times, 'y': np.full(3, 0.25)},
        'tha_psif_velocity_gsm': {
            'x': times,
            'y': np.tile([500.0, 0.0, 0.0], (3, 1)),
        },
        'tha_psif_avgtemp': {'x': times, 'y': np.full(3, 10000.0)},
    }

    moments = vis.time_averaged_particle_moments(
        esa, sst, prefix='tha'
    )

    assert moments['density_cm3'] == pytest.approx(1.25)
    np.testing.assert_allclose(
        moments['velocity_gsm_kms'], [180.0, 0.0, 0.0]
    )
    assert moments['temperature_eV'] > 2800.0
    assert moments['instrument_mixing_temperature_eV'] > 0.0
    assert set(moments['instruments']) == {'esa', 'sst'}


def test_spectrum_rejects_invalid_or_empty_input():
    esa, _ = _products()
    with pytest.raises(ValueError, match='no samples'):
        vis.spectrum(esa, 'peif_en_eflux', prefix='tha', trange=['2025-01-01', '2025-01-02'])
    esa['tha_peif_en_eflux']['y'][:] = 0
    with pytest.raises(ValueError, match='no positive'):
        vis.spectrum(esa, 'peif_en_eflux', prefix='tha')


def test_get_ma_beta_interpolates_across_invalid_source_samples():
    magnetic_time = pd.date_range('2024-01-01', periods=5, freq='min')
    plasma_time = pd.date_range('2024-01-01', periods=3, freq='2min')
    fgm = {'tha_fgs_btotal': {'x': magnetic_time, 'y': np.full(5, 10.0)}}
    esa = {
        'tha_peir_density': {'x': plasma_time, 'y': np.array([0.1, np.nan, 0.3])},
        'tha_peir_avgtemp': {'x': plasma_time, 'y': np.array([100.0, np.nan, 300.0])},
    }
    mom = {
        'tha_peim_velocity_gsm': {
            'x': plasma_time,
            'y': np.array([[100.0, 0.0, 0.0], [np.nan, np.nan, np.nan], [300.0, 0.0, 0.0]]),
        }
    }

    _, mach_alfven, beta = vis.get_ma_beta(
        fgm, esa, mom, prefix='tha', thres_ma=None, thres_bt=None
    )

    assert np.isfinite(mach_alfven).all()
    assert np.isfinite(beta).all()


def test_alfven_mach_number_converts_kilometres_per_second_once():
    magnetic_field_nt = np.array([10.0])
    density_cm3 = np.array([0.1])
    speed_kms = np.array([100.0])
    expected = speed_kms * 1000.0 / plasma.alfven_velocity(
        magnetic_field_nt * 1e-5, density_cm3, 1
    )

    np.testing.assert_allclose(
        vis.th_MA(magnetic_field_nt, density_cm3, speed_kms), expected
    )


def test_velocity_parallel_perp_interpolates_field_to_mom_time():
    magnetic_time = pd.date_range('2025-01-01', periods=3, freq='min')
    velocity_time = pd.to_datetime(['2025-01-01 00:00:30', '2025-01-01 00:01:30'])
    fgm = {
        'thb_fgs_gsm': {
            'x': magnetic_time,
            'y': np.array([[10.0, 0.0, 0.0]] * 3),
        }
    }
    mom = {
        'thb_peim_velocity_gsm': {
            'x': velocity_time,
            'y': np.array([[3.0, 4.0, 0.0], [-3.0, 0.0, 4.0]]),
        }
    }

    result_time, parallel, perpendicular = vis.velocity_parallel_perp(fgm, mom)

    assert result_time.equals(pd.DatetimeIndex(velocity_time))
    np.testing.assert_allclose(parallel, [3.0, -3.0])
    np.testing.assert_allclose(perpendicular, [4.0, 4.0])

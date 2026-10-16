"""MMS product panels, built on plots.vis plotting primitives."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from plots import vis
from .timeseries import tplot as _tplot


@vis.pass_axes(subplot_spec=(3, 1, True, False))
def fpi(fpi, f=None, ax=None):
    _tplot(fpi, 'numberdensity_fast', f=f, ax=ax[0])
    ax[0].set_yscale('log')
    ax[0].set_ylabel(r'$N_i (cm^{-3})$')

    _tplot(fpi, 'tempperp_fast', f=f, ax=ax[1])
    _tplot(fpi, 'temppara_fast', f=f, ax=ax[1])
    ax[1].set_yscale('log')
    ax[1].set_ylabel(r'$T (eV)$')

    _tplot(fpi, 'bulkv_gse_fast', f=f, ax=ax[2])
    ax[2].set_ylabel(r'$V_i (km/s)$')
    for i in range(2):
        ax[i + 1].legend()
    return f, ax


@vis.pass_axes(subplot_spec=(5, 1, True, False))
def hpca(hpca, f=None, ax=None):
    _tplot(hpca, 'hplus_number_density', f=f, ax=ax[0])
    _tplot(hpca, 'heplus_number_density', f=f, ax=ax[0])
    _tplot(hpca, 'oplus_number_density', f=f, ax=ax[0])
    ax[0].set_yscale('log')
    ax[0].set_ylabel(r'$N_i (cm^{-3})$')

    _tplot(hpca, 'hplus_tperp', f=f, ax=ax[1])
    _tplot(hpca, 'heplus_tperp', f=f, ax=ax[1])
    _tplot(hpca, 'oplus_tperp', f=f, ax=ax[1])
    ax[1].set_yscale('log')
    ax[1].set_ylabel(r'$T_\perp (eV)$')

    _tplot(hpca, 'hplus_ion_bulk_velocity_GSM', f=f, ax=ax[2])
    _tplot(hpca, 'heplus_ion_bulk_velocity_GSM', f=f, ax=ax[3])
    _tplot(hpca, 'oplus_ion_bulk_velocity_GSM', f=f, ax=ax[4])
    ax[2].set_ylabel(r'$V_i (km/s)$')
    ax[3].set_ylabel(r'$V_i (km/s)$')
    ax[4].set_ylabel(r'$V_i (km/s)$')
    for i in range(5):
        ax[i].legend()
    return f, ax

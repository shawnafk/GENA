from merge import gio
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta
from poliastro.twobody.sampling import EpochsArray
from scipy.interpolate import interp1d
from poliastro.core.perturbations import J2_perturbation
from poliastro.core.perturbations import third_body
from poliastro.core.propagation import func_twobody
from poliastro.twobody.propagation import CowellPropagator
from poliastro.ephem import build_ephem_interpolant
from astropy import units as u
from poliastro.bodies import Moon, Earth, Sun
from poliastro.twobody import Orbit




Moon_J2 = 2.036e-4
k_moon = Moon.k.to(u.km**3 / u.s**2).value
k_earth = Earth.k.to(u.km**3 / u.s**2).value
k_sun = Earth.k.to(u.km**3 / u.s**2).value

epoch_start = Time("2024-07-05 00:00", scale="utc")
epoch_end = Time("2025-12-05 00:00", scale="utc")
body_moon_2_earth = build_ephem_interpolant(
    Moon,
    28 * u.day,
    (epoch_start, epoch_end),
    rtol=1e-2,
)
body_sun_2_earth = build_ephem_interpolant(
    Sun,
    365 * u.day,
    (epoch_start, epoch_end),
    rtol=1e-2,
)
def f_pos_earth(t):
    return -1 * body_moon_2_earth(t)
def f_pos_sun(t):
    return body_sun_2_earth(t) - body_moon_2_earth(t)


def f(t0, u_, k):
    du_kep = func_twobody(t0, u_, k)
    a_moon_j2 = J2_perturbation(
        t0, u_, k, J2=Moon_J2, R=Moon.R.to(u.km).value
    )
    a_earth = third_body(
        t0, u_, k_moon, k_earth, f_pos_earth
    )
    a_sun = third_body(
        t0, u_, k_moon, k_sun, f_pos_sun
    )
    sum_a = a_moon_j2 + a_earth + a_sun
    du_ad = np.array([0, 0, 0, sum_a[0], sum_a[1], sum_a[2]])
    return du_kep + du_ad


pdf = gio.from_pqt('pdf.pqt')
t,r = gio.get_pos(pdf)
t_obs=t[::60]
r_obs=r[::60]/1000
v_obs = np.gradient(r_obs,axis=0)/np.gradient(t_obs,axis=0)[:,np.newaxis]


fx = interp1d(t_obs,r_obs[:,0])
fy = interp1d(t_obs,r_obs[:,1])
fz = interp1d(t_obs,r_obs[:,2])

fvx = interp1d(t_obs,v_obs[:,0])
fvy = interp1d(t_obs,v_obs[:,1])
fvz = interp1d(t_obs,v_obs[:,2])

t_str = "2024-07-10 20:50:00"
t_timestamp = pd.to_datetime(t_str).timestamp() 
r_init = [fx(t_timestamp),fy(t_timestamp),fz(t_timestamp)] << u.km
v_init = [fvx(t_timestamp),fvy(t_timestamp),fvz(t_timestamp)] << u.km / u.s


orb = Orbit.from_vectors(Moon, r_init, v_init,  epoch=Time(t_str,scale="utc"))


tofs = TimeDelta(np.linspace(0, 30 * u.day, num=10000))
ephem = orb.to_ephem(
    EpochsArray(orb.epoch + tofs, method=CowellPropagator(rtol=1e-6, f=f)),
)

r_predict,v_predict = ephem.rv()

fig,ax = plt.subplots(3,1,sharex=True)
for i in range(3):
    ax[i].plot(pd.to_datetime(t_obs[::100],unit='s'),r_obs[::100,i],'x',label='measured')
    ax[i].plot(pd.to_datetime(ephem.epochs.value),r_predict[:,i],label='predicted')
    ax[i].set_ylabel(['x (km)', 'y (km)', 'z (km)'][i])
    ax[i].legend()
ax[2].set_xlabel('Time')
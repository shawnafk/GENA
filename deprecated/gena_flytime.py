import numpy as np

def tof(E,mass=1,l=50e-3):
    mH = 1836 * 511000
    c = 2.997e8
    v = np.sqrt(2*E/(mass*mH))*c
    t = l/v * 1e9
    return t

#H
E = np.logspace(-1,3)*1000
tof_H = tof(E)
tof_O = tof(E,15.872162485065711)


from myplot.fig_head import *
import matplotlib.pyplot as plt

f,ax = plt.subplots(1,1,figsize=(6,4))
ax.loglog(E/1000,tof_H, label='H', color='red')
ax.loglog(E/1000,tof_O, label='O', color='blue')

ax.axhspan(1, 25, color='red', alpha=0.3)
ax.axhspan(114, 1000, color='blue', alpha=0.3)
ax.legend()
ax.set_xlabel('E (keV)')
ax.set_ylabel('ToF (ns)')
ax.set_xlim(1,300)
ax.set_ylim(1,1000)
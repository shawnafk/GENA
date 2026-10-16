import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#use SI unit
G = 6.67408e-11 #m^3 kg^-1 s^-2
M_earth = 5.972e24
M_H = 1.6735575e-27
RE = 6371000
q_H = 1.60217662e-19
kB = 1.38064852e-23
speed_of_light = 299792458  #m/s
#y is coordinates and momentum [r, p]
#Hamiltonian is H(r,p) = p^2/2m - GMm/r
#dr/dt = p/m
#dp/dt = - Gmm/r^2
def gravity(t, var):
    x, y, px, py = var
    r = np.sqrt(x**2 + y**2)
    dxdt = px / M_H
    dydt = py / M_H
    dpxdt = -G * M_earth * M_H / r**2  * x/r
    dpydt = -G * M_earth * M_H / r**2  * y/r
    dydt = [dxdt, dydt, dpxdt, dpydt]
    return dydt


#assuming E in 10 - 100 eV
#E = np.linspace(10,100,1) * q_H * 
#p = np.sqrt(2*M_H*E)
v0 =10000*M_H
alpha = np.pi/4
y0 = [60*RE, 0, v0*cos(alpha), v0*sin(alpha) ]

t_span = [0, 20000]

sol = solve_ivp(gravity, t_span, y0, t_eval = np.arange(0, t_span[-1], 1))
H = (sol.y[2]**2 + sol.y[3]**2)/2/M_H - G*M_earth*M_H/np.sqrt(sol.y[0]**2 + sol.y[1]**2)

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from math import pi, sqrt

radius = 6378.0 # km
mu = 3.986E+05 # km^3/s^2

# y is a vector: [rx ry rz vx vy vz]
def f(t, y):
    rx = y[0]
    ry = y[1]
    rz = y[2]
    vx = y[3]
    vy = y[4]
    vz = y[5]

    d = (rx**2+ry**2+rz**2)**(3/2)

    dvx_dt = -mu*rx/d
    dvy_dt = -mu*ry/d
    dvz_dt = -mu*rz/d
    drx_dt = vx
    dry_dt = vy
    drz_dt = vz

    return np.array([drx_dt, dry_dt, drz_dt, dvx_dt, dvy_dt, dvz_dt])

#r0_norm = radius + 600.0 # km
T_new = 24*60*60 # 1 day in seconds
R_new = ((mu*T_new**2)/(4*pi**2))**(1/3)

print("R_new=",R_new)


#T = 2*pi*sqrt((r0_norm**3)/mu) # seconds
c = 10
times_per_T = 30
t_span = np.array([0, c*T_new])
times = np.linspace(t_span[0], t_span[1], times_per_T*c)

v0_norm = sqrt(mu/R_new) # velocity for circular motion

print(v0_norm)

y0 = np.array([R_new, 0, 0, 0, v0_norm, 0])

sol = solve_ivp(f, t_span, y0, t_eval=times, method='LSODA')
t = sol.t
rx = sol.y[0]
ry = sol.y[1]
rz = sol.y[2]
vx = sol.y[3]
vy = sol.y[4]
vz = sol.y[5]

#print(T)

# plot of answer
fig, ax = plt.subplots(7,1,sharex=True)
ax[0].plot(t, vx, 'r-', label='velocity x')
ax[0].set_ylabel('vx (m/s)')
ax[0].legend()
ax[1].plot(t, vy, 'g-', label='velocity y')
ax[1].set_ylabel('vy (m/s)')
ax[1].legend()
ax[2].plot(t, vz, 'b-', label='velocity z')
ax[2].set_ylabel('vz (m/s)')
ax[2].legend()

ax[3].plot(t, rx, 'g-', label='position x')
ax[3].plot([t[0], t[-1]], [0,0], 'k--')
ax[3].set_ylabel('rx (m)')
ax[3].legend()
ax[4].plot(t, ry, 'g-', label='position y')
ax[4].plot([t[0], t[-1]], [0,0], 'k--')
ax[4].set_ylabel('ry (m)')
ax[4].legend()
ax[5].plot(t, rz, 'g-', label='position z')
ax[5].plot([t[0], t[-1]], [0,0], 'k--')
ax[5].set_ylabel('rz (m)')
ax[5].legend()

ax[6].plot(rx, ry)
plt.show()

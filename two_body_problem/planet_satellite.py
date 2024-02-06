import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

radius = 6378.0 # km
mu = 3.986E+05 # km^3/s^2

def f(t, y):
    vx = y[0]
    vy = y[1]
    vz = y[2]
    rx = y[3]
    ry = y[4]
    rz = y[5]

    d = (rx**2+ry**2+rz**2)**(3/2)

    dvx_dt = -mu*rx/d
    dvy_dt = -mu*ry/d
    dvz_dt = -mu*rz/d
    drx_dt = vx
    dry_dt = vy
    drz_dt = vz

    return np.array([dvx_dt, dvy_dt, dvz_dt, drx_dt, dry_dt, drz_dt])

t_span = np.array([0, 100*60])
#times = np.linspace(t_span[0], t_span[1], 101)

r0_norm = radius + 450.0 # km
v0_norm = (mu/radius)**0.5
y0 = np.array([0, v0_norm, 0, r0_norm, 0, 0])

#sol = solve_ivp(f, t_span, y0, t_eval=times)
sol = solve_ivp(f, t_span, y0)
t = sol.t
vx = sol.y[0]
vy = sol.y[1]
vz = sol.y[2]
rx = sol.y[3]
ry = sol.y[4]
rz = sol.y[5]

# plot of answer
fig, ax = plt.subplots(6,1,sharex=True)
ax[0].plot(t, vx, 'r-', label='velocity x')
ax[0].set_ylabel('vx (m/s)');
ax[0].legend()
ax[1].plot(t, vy, 'g-', label='velocity y')
ax[1].set_ylabel('vy (m/s)');
ax[1].legend()
ax[2].plot(t, vz, 'b-', label='velocity z')
ax[2].set_ylabel('vz (m/s)');
ax[2].legend()

ax[3].plot(t, rx, 'g-', label='position x')
ax[3].plot([t[0], t[-1]], [0,0], 'k--')
ax[3].set_ylabel('rx (m)');
ax[3].legend()
ax[4].plot(t, ry, 'g-', label='position y')
ax[4].plot([t[0], t[-1]], [0,0], 'k--')
ax[4].set_ylabel('ry (m)');
ax[4].legend()
ax[5].plot(t, rz, 'g-', label='position z')
ax[5].plot([t[0], t[-1]], [0,0], 'k--')
ax[5].set_ylabel('rz (m)')
ax[5].legend()
plt.show()

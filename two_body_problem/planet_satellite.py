import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from math import pi, sqrt

# Earth values
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

def get_solution(t_span, y0, times, method='Radau'):
    return solve_ivp(f, t_span, y0, t_eval=times, method=method)

def graph_solution(solution):
    t = solution.t
    y = []
    for i in range(6):
        y.append(solution.y[i])

    fig, ax = plt.subplots(6,1,sharex=True)
    ax[0].plot(t, y[3], 'r-', label='velocity x')
    ax[0].set_ylabel('vx (m/s)')
    ax[0].legend()
    ax[1].plot(t, y[4], 'g-', label='velocity y')
    ax[1].set_ylabel('vy (m/s)')
    ax[1].legend()
    ax[2].plot(t, y[5], 'b-', label='velocity z')
    ax[2].set_ylabel('vz (m/s)')
    ax[2].legend()

    ax[3].plot(t, y[0], 'g-', label='position x')
    ax[3].plot([t[0], t[-1]], [0,0], 'k--')
    ax[3].set_ylabel('rx (m)')
    ax[3].legend()
    ax[4].plot(t, y[1], 'g-', label='position y')
    ax[4].plot([t[0], t[-1]], [0,0], 'k--')
    ax[4].set_ylabel('ry (m)')
    ax[4].legend()
    ax[5].plot(t, y[2], 'g-', label='position z')
    ax[5].plot([t[0], t[-1]], [0,0], 'k--')
    ax[5].set_ylabel('rz (m)')
    ax[5].legend()

    fig2, ax2 = plt.subplots()
    ax2.plot(y[0], y[1])

    plt.show()

def get_cyclic_initial_values(c=10):
    #r0_norm = radius + 600.0 # km
    T_new = 24*60*60 # 1 day in seconds
    R_new = ((mu*T_new**2)/(4*pi**2))**(1/3)

    print("T_new =",T_new)
    print("R_new =",R_new)

    #T = 2*pi*sqrt((r0_norm**3)/mu) # seconds
    #c = 100
    times_per_T = 30
    t_span = np.array([0, c*T_new])
    times = np.linspace(t_span[0], t_span[1], times_per_T*c)

    v0_norm = sqrt(mu/R_new) # velocity for circular motion

    print("v0_norm =",v0_norm)

    y0 = np.array([R_new, 0, 0, 0, v0_norm, 0])

    # test calculate orbital period
    energy = 0.5*v0_norm**2 - mu/R_new
    a = -mu/(2*energy)
    orbital_T = 2*pi*sqrt((a**3)/mu)
    print("orbital_T =",orbital_T)

    return t_span, times, y0

def get_general_initial_values(r0, v0, c=10, times_per_T=30):
    r0_norm = np.linalg.norm(r0)
    v0_norm = np.linalg.norm(v0)

    energy = 0.5*v0_norm**2 - mu/r0_norm
    a = -mu/(2*energy)
    print("a =",a)
    T = 2*pi*sqrt((a**3)/mu)

    t_span = np.array([0, c*T])
    times = np.linspace(t_span[0], t_span[1], times_per_T*c)
    y0 = np.array([r0[0], r0[1], r0[2], v0[0], v0[1], v0[2]])
    
    return t_span, times, y0

def test_solution(graph=True):
    #t_span, times, y0 = get_cyclic_initial_values(c=1)

    h = 600 # km
    r0 = np.array([radius+h, 0, 0])

    v_cyclic = sqrt(mu/r0[0]) # velocity for circular motion
    print("v_cyclic =",v_cyclic)
    v_escape = sqrt(2*mu/r0[0]) # escape velocity
    print("v_escape =",v_escape)

    v0 = np.array([0, v_escape-0.001, 0])

    t_span, times, y0 = get_general_initial_values(r0=r0, v0=v0, c=5, times_per_T=50)
    sol = get_solution(t_span=t_span, y0=y0, times=times)
    if(graph): 
        graph_solution(sol)

#test_solution()
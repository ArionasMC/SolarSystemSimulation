import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from math import pi, sqrt

m1 = 5.972e24 # earth (kg)
m2 = m1 # second earth
G = 6.6743e-20 # Gravitational Constant (km^3 * kg^-1 * s^-2)
M = m1+m2
mu = G*M # (km^3 * s^-2)
mu1 = G*m1
mu2 = G*m2
R = 6378.0 # km

# state vector y = [r1x r1y r1z v1x v1y v1z r2x r2y r2z v2x v2y v2z]
def f(t, y):
    y1 = y[:6]
    y2 = y[6:]
    r = y2[:3]-y1[:3]
    r_norm = np.linalg.norm(r)

    dv1_dt = (mu2/r_norm**3) * r
    dv2_dt = -(mu1/r_norm**3) * r
    dr1_dt = y1[3:6] # v1 vector
    dr2_dt = y2[3:6] # v2 vector
    
    result = np.concatenate((dr1_dt, dv1_dt, dr2_dt, dv2_dt))
    return result

PLOT_NAMES = ('position x', 'position y', 'position z',
              'velocity x', 'velocity y', 'velocity z')
LABEL_NAMES = ('$r_x$ (m)', '$r_y$ (m)', '$r_z$ (m)',
               '$v_x$ (m/s)', '$v_y$ (m/s)', '$v_z$ (m/s)')

def graph_solution(solution):
    t = solution.t
    y = solution.y

    for i in range(2):
        fig, ax = plt.subplots(6, 1, sharex=True)
        for j in range(6):
            ax[j].plot(t, y[i*6+j], 'r', label=PLOT_NAMES[j])
            ax[j].set_ylabel(LABEL_NAMES[j])
            ax[j].legend()
        ax[5].set_xlabel('Time (s)')
        fig.suptitle('Body m'+str(i+1))
    
    fig, ax = plt.subplots()
    ax.plot(y[0], y[1], 'b', label='$r_{1y}$($r_{1x}$)')
    ax.plot(y[6], y[7], 'r', label='$r_{2y}$($r_{2x}$)')
    ax.set_xlabel('Position X')
    ax.set_ylabel('Position Y')
    ax.legend()

    plt.show()

# executes only if the file is run
if __name__ == '__main__':
    t_span = np.array([0, 60*60*24]) # 1 day in seconds
    times = np.linspace(t_span[0], t_span[1], 101)

    # circular initial values
    # let R1 = (0, 0, 0) so R2 = r at the beginning
    r0 = R + 600 # km
    v0 = sqrt(mu/r0)
    state1 = np.array([0, 0, 0, 0, 0, 0])
    state2 = np.array([r0, 0, 0, 0, v0, 0])
    y0 = np.append(state1, state2)

    solution = solve_ivp(f, t_span, y0, t_eval=times, method='Radau')
    graph_solution(solution)
from ursina import *
from ursina.shaders import lit_with_shadows_shader
from ursina.prefabs.first_person_controller import FirstPersonController
from core.celestial_body import CelestialBody
from two_body_problem.planet_satellite import *
import numpy as np

app = Ursina()

camera = FirstPersonController(gravity=0, z=-3, y=-2)

earth = CelestialBody("Earth", mass=100.0, model="sphere", color=color.blue, scale=1, shader=lit_with_shadows_shader)
satellite = CelestialBody("Satellite", mass=0.1, model="sphere", color=color.gray, scale=0.1, shader=lit_with_shadows_shader)

TEST_SCALE = radius
h = 600 #km
r0 = np.array([radius+h, 0, 0])
satellite.position = r0/TEST_SCALE
print(r0/TEST_SCALE)

v_cyclic = sqrt(mu/r0[0]) # velocity for circular motion
v_escape = sqrt(2*mu/r0[0]) # escape velocity
v0 = np.array([0, v_escape*0.8, 0])

times_per_T = 50
c = 10
t_span, times, y0 = get_general_initial_values(r0=r0, v0=v0, c=c, times_per_T=times_per_T)
sol = get_solution(t_span=t_span, y0=y0, times=times)
#graph_solution(sol)
print("sol len =",len(sol.y[0]),len(sol.y[1]),len(sol.y[2]))

index = 0

speed_y = 2
bodies = [earth, satellite]

for body in bodies:
    print(body.name, body.position)

def update():
    global index
    camera.y += held_keys['space'] * time.dt * speed_y
    camera.y -= held_keys['shift'] * time.dt * speed_y

    satellite.position = np.array([sol.y[0][index], sol.y[2][index], sol.y[1][index]])/TEST_SCALE
    index += 1
    index %= c*times_per_T

def input(key):
    #print(key)
    if key == 'escape':
        camera.enabled = not(camera.enabled)

app.run()
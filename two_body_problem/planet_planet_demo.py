from ursina import *
from ursina.shaders import lit_with_shadows_shader
from ursina.prefabs.first_person_controller import FirstPersonController
from .core import celestial_body
from planet_planet import *
import numpy as np

app = Ursina()

camera = FirstPersonController(gravity=0, z=-3, y=-2)

earth = celestial_body.CelestialBody("Earth", mass=100.0, model="sphere", color=color.blue, scale=1, shader=lit_with_shadows_shader)
earth2 = celestial_body.CelestialBody("Earth2", mass=100.0, model="sphere", color=color.red, scale=1, shader=lit_with_shadows_shader)

TEST_SCALE = R
h = 600 #km
r0 = np.array([R+h, 0, 0])
earth2.position = r0/TEST_SCALE
print(r0/TEST_SCALE)

v_cyclic = sqrt(mu/r0[0]) # velocity for circular motion
v0 = np.array([0, v_cyclic, 0])

t_span = np.array([0, 60*60*24]) # 1 day in seconds
times_per_T = 101
c = 1
times = np.linspace(t_span[0], t_span[1], times_per_T)
state1 = np.array([0, 0, 0, 0, 0, 0])
state2 = np.array([r0[0], 0, 0, 0, v0[1], 0])
y0 = np.append(state1, state2)

sol = solve_ivp(f, t_span, y0, t_eval=times, method='Radau')

index = 0

speed_y = 2
bodies = [earth, earth2]

for body in bodies:
    print(body.name, body.position)

def update():
    global index
    camera.y += held_keys['space'] * time.dt * speed_y
    camera.y -= held_keys['shift'] * time.dt * speed_y

    earth2.position = np.array([sol.y[0][index], sol.y[2][index], sol.y[1][index]])/TEST_SCALE
    index += 1
    index %= c*times_per_T

def input(key):
    #print(key)
    if key == 'escape':
        camera.enabled = not(camera.enabled)

app.run()
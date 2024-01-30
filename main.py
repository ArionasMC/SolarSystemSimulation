from ursina import *
from ursina.shaders import lit_with_shadows_shader
from ursina.prefabs.first_person_controller import FirstPersonController
from core.celestial_body import CelestialBody

app = Ursina()

camera = FirstPersonController(gravity=0, z=-3, y=-2)

sun = CelestialBody(mass=1.0, model="sphere", color=color.yellow, scale=2, shader=lit_with_shadows_shader)
light = PointLight(parent=sun, color=color.white)
light.position = (0, 0, 0)

planet = CelestialBody(mass=0.1, model="sphere", color=color.blue, scale=0.5, shader=lit_with_shadows_shader)
planet.position = (-2,0,0)

speed_y = 2

def update():
    camera.y += held_keys['space'] * time.dt * speed_y
    camera.y -= held_keys['shift'] * time.dt * speed_y

def input(key):
    #print(key)
    if key == 'escape':
        camera.enabled = not(camera.enabled)

app.run()
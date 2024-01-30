from ursina import *
from ursina.shaders import lit_with_shadows_shader
from ursina.prefabs.first_person_controller import FirstPersonController
from core.celestial_body import CelestialBody

app = Ursina()

camera = FirstPersonController(gravity=0, z=-3, y=-2)

sun = CelestialBody("Sun", mass=100.0, model="sphere", color=color.yellow, scale=2, shader=lit_with_shadows_shader)
light = PointLight(parent=sun, color=color.white)
light.position = (0, 0, 0)

planet = CelestialBody("Planet 1", mass=0.1, model="sphere", color=color.blue, scale=0.5, shader=lit_with_shadows_shader)
planet.position = (-2,0,0)
planet.velocity = [0,0,0] #it works, must set otherwise sun moves too (?!!??!)
planet.acceleration = [0,0,0]

speed_y = 2
bodies = [sun, planet]

for body in bodies:
    print(body.name, body.acceleration)

def calculate_gravity(body1, body2):
    dx = body2.x - body1.x
    dy = body2.y - body1.y
    distance = math.sqrt(dx**2 + dy**2)
    force = 6.67430e-11 * (body1.mass * body2.mass) / distance**2
    angle = math.atan2(dy, dx)
    fx = force * math.cos(angle)
    fy = force * math.sin(angle)
    return fx, fy



def update():
    camera.y += held_keys['space'] * time.dt * speed_y
    camera.y -= held_keys['shift'] * time.dt * speed_y

    for body in bodies:
        body.update_vectors()

    fx, fz = calculate_gravity(sun, planet)
    planet.acceleration = (fx, 0, fz)

def input(key):
    #print(key)
    if key == 'escape':
        camera.enabled = not(camera.enabled)

app.run()
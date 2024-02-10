from ursina import Entity, time

class CelestialBody(Entity):
    def __init__(self, name, mass=1.0, velocity=[0, 0, 0], acceleration=[0,0,0], add_to_scene_entities=True, **kwargs):
        super().__init__(add_to_scene_entities, **kwargs)
        self.name = name
        self.mass = mass
        self.velocity = velocity # (vx, vy, vz)
        self.acceleration = acceleration # (ax, ay, az)

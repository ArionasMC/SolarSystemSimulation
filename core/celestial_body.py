from ursina import Entity, time

class CelestialBody(Entity):
    def __init__(self, name, mass=1.0, velocity=[0, 0, 0], acceleration=[0,0,0], add_to_scene_entities=True, **kwargs):
        super().__init__(add_to_scene_entities, **kwargs)
        self.name = name
        self.mass = mass
        self.velocity = velocity # (vx, vy, vz)
        self.acceleration = acceleration # (ax, ay, az)

    def update_vectors(self):
        for i in range(3):
            self.velocity[i] += self.acceleration[i] * time.dt
        self.x += self.velocity[0] * time.dt
        self.y += self.velocity[1] * time.dt
        self.z += self.velocity[2] * time.dt
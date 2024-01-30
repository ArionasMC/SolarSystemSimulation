from ursina import Entity

class CelestialBody(Entity):
    def __init__(self, mass=1.0, add_to_scene_entities=True, **kwargs):
        super().__init__(add_to_scene_entities, **kwargs)
        self.mass = mass
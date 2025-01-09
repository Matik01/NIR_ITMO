class Protein:
    def __init__(self, coords, kdtrees):
        self._coords = coords
        self._kdtrees = kdtrees
    def get_coords(self):
        return self._coords
    def get_kdtrees(self):
        return self._kdtrees
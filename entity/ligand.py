import numpy as np


class Ligand:
    def __init__(self, atom_types: list, coords: np.array):
        self._atom_types = atom_types
        self._coords = coords

    def get_coords(self):
        return self._coords

    def get_atom_types(self):
        return self._atom_types

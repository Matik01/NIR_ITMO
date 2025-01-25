import numpy as np


class Ligand:
    def __init__(self, atom_types: list, coords: np.array):
        self._atom_types: list = atom_types
        self._coords: np.array = coords

    def get_coords(self):
        return self._coords

    def get_atom_types(self):
        return self._atom_types

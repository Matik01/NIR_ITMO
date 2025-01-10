import numpy as np
from sklearn.neighbors import KDTree


class Protein:
    def __init__(self, atom_types: list, coords: dict, kd_trees: dict):
        self._atom_types = atom_types
        self._coords = coords
        self._kd_trees = kd_trees

    def get_coords(self):
        return self._coords

    def get_kd_trees(self):
        return self._kd_trees

    def get_atom_types(self):
        return self._atom_types

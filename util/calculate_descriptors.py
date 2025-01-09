from entity.ligand import Ligand
from entity.protein import Protein
from entity.constants import Constants
from sklearn.neighbors import KDTree
import numpy as np
from itertools import product


def read_pdb(filename: str):
    pdb_data = None
    with open(filename, 'r') as f_in:
        print(f'Reading file: {filename}')
        pdb_data = f_in.readlines()
    protein_data = [line.rstrip('\n') for line in pdb_data if line[17:20] in Constants.AMINOACID_LIST]
    ligand_data = [line.rstrip('\n') for line in pdb_data if line[17:20] == Constants.LIGAND]



def fc_cutoff(r_distance: float) -> float:
    if r_distance <= Constants.RC_CUTOFF_DISTANCE:
        return np.tanh(1 - r_distance/Constants.RC_CUTOFF_DISTANCE) ** 3
    else:
        return 0

def distance_function(eta: float, r_distance: float, r_s: float) -> float:
    return np.sum(np.exp(-eta*(r_distance - r_s)**2) * fc_cutoff(r_distance)) # i != j

def angle_function(zeta: float, eta: float,
                   lambda_values: tuple[1, -1], cos_teta_ijk: float,
                   distance_ij: float, distance_jk: float, distance_ik: float) -> list:
    angle_descriptors = []
    for lambda_value in lambda_values:
        g_angular = 2**(1 - zeta) * np.sum((1 + lambda_value * cos_teta_ijk)**zeta * np.exp(-eta*(distance_ij**2 + distance_jk**2 + distance_ik**2)) * fc_cutoff(distance_ij) * fc_cutoff(distance_jk) * fc_cutoff(distance_ik))
        angle_descriptors.append(g_angular)
    return angle_descriptors

if __name__ == '__main__':
    read_pdb('test_system.pdb')

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
    mol = process_ligand(ligand_data)
    process_protein(protein_data)


def process_ligand(ligand_data: list):
    coords = np.array([]).reshape((0, 3)).astype(float)
    atom_types = []
    '''
    Getting atom types and coords for each
    '''
    for item in ligand_data:
        atom_type = item[76:78]
        atom_types.append(atom_type.strip())
        x = item[32:38]
        y = item[40:46]
        z = item[48:54]
        coords = np.vstack([coords, [x, y, z]]).astype(float)
    return Ligand(atom_types=atom_types, coords=coords)


def process_protein(protein_data: list):
    coords = {atom_type: np.array([]).reshape((0,3)).astype(float) for atom_type in Constants.STANDARD_ATOMS_MAP.keys()}
    atom_types = []
    for item in protein_data:
        atom_type = item[76:78].strip()
        atom_types.append(atom_type)
        x = item[32:38]
        y = item[40:46]
        z = item[48:54]
        coords[atom_type] = np.vstack([coords[atom_type], [x, y, z]]).astype(float)

    kd_trees = dict()
    for atom_type in atom_types:
        if len(coords[atom_type]) != 0:
            kd_trees[atom_type] = KDTree(coords[atom_type])

    return Protein(atom_types=atom_types, coords=coords, kd_trees=kd_trees)



def fc_cutoff(r_distance: float) -> float:
    if r_distance <= Constants.RC_CUTOFF_DISTANCE:
        return np.tanh(1 - r_distance / Constants.RC_CUTOFF_DISTANCE) ** 3
    else:
        return 0


def distance_function(eta: float, r_distance: float, r_s: float) -> float:
    return np.sum(np.exp(-eta * (r_distance - r_s) ** 2) * fc_cutoff(r_distance))  # i != j


def angle_function(zeta: float, eta: float,
                   lambda_values: tuple[1, -1], cos_teta_ijk: float,
                   distance_ij: float, distance_jk: float, distance_ik: float) -> list:
    angle_descriptors = []
    for lambda_value in lambda_values:
        g_angular = 2 ** (1 - zeta) * np.sum((1 + lambda_value * cos_teta_ijk) ** zeta * np.exp(
            -eta * (distance_ij ** 2 + distance_jk ** 2 + distance_ik ** 2)) * fc_cutoff(distance_ij) * fc_cutoff(
            distance_jk) * fc_cutoff(distance_ik))
        angle_descriptors.append(g_angular)
    return angle_descriptors


if __name__ == '__main__':
    read_pdb('test_system.pdb')

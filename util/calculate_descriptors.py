from entity.ligand import Ligand
from entity.protein import Protein
from entity.constants import Constants
from sklearn.neighbors import KDTree
import numpy as np
from itertools import product
import pickle as pkl


def read_pdb(filename: str) -> (Ligand, Protein):
    pdb_data = None
    with open(filename, 'r') as f_in:
        print(f'Reading file: {filename}')
        pdb_data = f_in.readlines()
    protein_data = [line.rstrip('\n') for line in pdb_data if line[17:20] in Constants.AMINOACID_LIST]
    ligand_data = [line.rstrip('\n') for line in pdb_data if line[17:20] == Constants.LIGAND]

    return process_ligand(ligand_data), process_protein(protein_data)


def process_ligand(ligand_data: list) -> Ligand:
    coords = np.array([]).reshape((0, 3)).astype(float)
    atom_types = []
    '''
    Getting atom types and coords for each
    '''
    for item in ligand_data:
        atom_char = item[76:78]
        atom_types.append(atom_char.strip())
        x = item[32:38]
        y = item[40:46]
        z = item[48:54]
        coords = np.vstack([coords, [x, y, z]]).astype(float)
    return Ligand(atom_types=atom_types, coords=coords)


def process_protein(protein_data: list) -> Protein:
    coords = {atom_char: np.array([]).reshape((0, 3)).astype(float) for atom_char in
              Constants.STANDARD_ATOMS_MAP.keys()}
    atom_types = list()
    for item in protein_data:
        atom_char = item[76:78].strip()
        if atom_char not in atom_types:
            atom_types.append(atom_char)
        x = item[32:38]
        y = item[40:46]
        z = item[48:54]
        coords[atom_char] = np.vstack([coords[atom_char], [x, y, z]]).astype(float)

    kd_trees = dict()
    for atom_char in atom_types:
        if len(coords[atom_char]) != 0:
            kd_trees[atom_char] = KDTree(coords[atom_char])

    keys_to_delete = set()
    for key in coords.keys():
        if coords[key].size == 0:
            keys_to_delete.add(key)

    if len(keys_to_delete) > 0:
        for key in keys_to_delete:
            print(f'Deleting atom char: {key}')
            del coords[key]

    return Protein(atom_types=atom_types, coords=coords, kd_trees=kd_trees)


def fc_cutoff(r_distance: float) -> float:
    return np.tanh(1 - r_distance / Constants.RC_CUTOFF_DISTANCE) ** 3


def distance_function(eta: float, distance_ij: float, r_s: float) -> float:
    return np.sum(np.exp(-eta * (distance_ij - r_s) ** 2) * fc_cutoff(distance_ij))  # i != j


def angle_function(zeta: float, eta: float, cosine_teta: float,
                   distance_ij: float, distance_jk: float, distance_ik: float, lambda_value: float) -> list:
    g_angular = 2 ** (1 - zeta) * np.sum((1 + lambda_value * cosine_teta) ** zeta * np.exp(
        -eta * (distance_ij ** 2 + distance_jk ** 2 + distance_ik ** 2)) * fc_cutoff(distance_ij) * fc_cutoff(
        distance_jk) * fc_cutoff(distance_ik))

    return g_angular


def bps_compute(prot: Protein, ligand: Ligand) -> list:
    env_descriptors = list()
    '''
    Output structure by index:
    0 - protein atom
    1 - ligand atom
    2 - calculated descriptors
    '''
    for atom_char in prot.get_coords().keys():
        idx, distances = prot.get_kd_trees()[atom_char].query_radius(ligand.get_coords(),
                                                                     Constants.RC_CUTOFF_DISTANCE,
                                                                     return_distance=True)
        descriptors_by_atom_id = list()
        for atom_id in range(ligand.get_coords().shape[0]):
            descriptors = list()
            env_coords = prot.get_coords()[atom_char][idx[atom_id], :]
            ligand_atom_coords = ligand.get_coords()[atom_id]

            j = np.tile(env_coords, (env_coords.shape[0], 1))
            k = np.repeat(env_coords, env_coords.shape[0], axis=0)
            distance_ij = np.sum((ligand_atom_coords - j) ** 2, axis=1) ** 0.5
            distance_ik = np.sum((ligand_atom_coords - k) ** 2, axis=1) ** 0.5
            distance_jk = np.sum((j - k) ** 2, axis=1) ** 0.5
            cosine_teta = np.sum((ligand_atom_coords - j) * (ligand_atom_coords - k), axis=1) / (
                    distance_ij * distance_ik)

            for result in cosine_teta:
                if int(result) > 1:
                    print(f'Error not valid cosine: {result}')

            for eta, rs in product(Constants.ETA_VALUES, Constants.RS_VALUES):
                temp = distance_function(eta=eta, distance_ij=distances[atom_id], r_s=rs)

                descriptors.append(temp)

            for lambda_value in Constants.LAMBDA_VALUES:
                for eta, zeta in product(Constants.ETA_VALUES, Constants.ZETA_VALUES):
                    temp = angle_function(zeta=zeta, eta=eta, cosine_teta=cosine_teta, distance_ij=distance_ij,
                                          distance_ik=distance_ik, distance_jk=distance_jk, lambda_value=lambda_value)
                    descriptors.append(temp)
            descriptors_by_atom_id.append(descriptors)
        env_descriptors.append(descriptors_by_atom_id)
    return env_descriptors


if __name__ == '__main__':
    files = ['test_system.pdb']  # TO ADD: Directory reader
    ligand_list: list[Ligand] = []
    protein_list: list[Protein] = []

    for file in files:  # Do I really need this?
        mol, prot = read_pdb(file)
        ligand_list.append(mol)
        protein_list.append(prot)

    for prot, ligand in zip(protein_list, ligand_list):
        out = bps_compute(prot=prot, ligand=ligand)
        out = np.array(out)
        np.savetxt('out.txt', out[0], delimiter=',')

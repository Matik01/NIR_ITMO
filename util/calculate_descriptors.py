from entity.ligand import Ligand
from entity.protein import Protein
from entity.constants import Constants
from sklearn.neighbors import KDTree
import numpy as np
from itertools import product
import logging

logger = logging.getLogger(__name__)


class BPS:
    def __init__(self, file: str):
        self._file: str = file

    def main(self) -> tuple:
        logging.basicConfig(filename="descriptors.log", level=logging.ERROR)

        ligand: Ligand
        protein: Protein

        ligand, protein = self._read_pdb(self._file)
        bps_out = self._bps_compute(prot=protein, ligand=ligand)
        bps_out = np.transpose(np.array(bps_out), (1, 0, 2))

        ligand_one_hot_encode = np.zeros((len(ligand.get_coords()), len(Constants.LIGAND_ATOM_TYPES)))
        for i in range(ligand.get_coords().shape[0]):
            ligand_one_hot_encode[i, Constants.LIGAND_ATOM_TYPES.index(ligand.get_atom_types()[i])] += 1.

        bps_out = bps_out.reshape(bps_out.shape[0], bps_out.shape[1] * bps_out.shape[2])
        bps_out = np.concatenate((ligand_one_hot_encode, bps_out), axis=1).astype(np.float32)
        return self._file, bps_out

    def _read_pdb(self, filename: str) -> (Ligand, Protein):
        with open(filename, 'r') as f_in:
            logger.info(f'Reading file: {filename}')
            pdb_data = f_in.readlines()

        protein_data = [line.rstrip('\n') for line in pdb_data if line[17:20] in Constants.AMINOACID_LIST]
        ligand_data = [line.rstrip('\n') for line in pdb_data if line[17:20] == Constants.LIGAND_IDENTIFIER]

        return self._process_ligand(ligand_data), self._process_protein(protein_data)

    def _process_ligand(self, ligand_data: list) -> Ligand:
        coords = np.array([]).reshape((0, 3)).astype(float)
        atom_types = []
        '''
        Getting atom types and coords for each
        '''
        for item in ligand_data:
            atom_char = item[76:78]
            atom_types.append(atom_char.upper().strip())
            x = item[32:38]
            y = item[40:46]
            z = item[48:54]
            coords = np.vstack([coords, [x, y, z]]).astype(float)
        return Ligand(atom_types=atom_types, coords=coords)

    def _process_protein(self, protein_data: list) -> Protein:
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
                logger.info(f'Deleting atom char: {key}')
                del coords[key]

        return Protein(atom_types=atom_types, coords=coords, kd_trees=kd_trees)

    def _fc_cutoff(self, r_distance: float) -> float:
        return np.tanh(1 - r_distance / Constants.RC_CUTOFF_DISTANCE) ** 3

    def _radial_distribution_function(self, eta: float, distance_ij: float, r_s: float) -> float:
        return np.sum(np.exp(-eta * (distance_ij - r_s) ** 2) * self._fc_cutoff(distance_ij))  # i != j

    def _angle_function(self, zeta: float, eta: float, cosine_theta: float,
                        distance_ij: float, distance_jk: float, distance_ik: float, lambda_value: float) -> list:
        g_angular = 2 ** (1 - zeta) * np.sum((1 + lambda_value * cosine_theta) ** zeta * np.exp(
            -eta * (distance_ij ** 2 + distance_jk ** 2 + distance_ik ** 2)) * self._fc_cutoff(
            distance_ij) * self._fc_cutoff(distance_jk) * self._fc_cutoff(distance_ik))

        return g_angular

    def _bps_compute(self, prot: Protein, ligand: Ligand) -> list:
        env_descriptors = list()
        '''
        Output array structure by index:
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
                cosine_theta = np.sum((ligand_atom_coords - j) * (ligand_atom_coords - k), axis=1) / (
                        distance_ij * distance_ik)

                for result in cosine_theta:
                    if int(result) > 1:
                        logger.error(f'Error not valid cosine: {result}')
                        return list()

                for eta, rs in product(Constants.ETA_VALUES, Constants.RS_VALUES):
                    temp = self._radial_distribution_function(eta=eta, distance_ij=distances[atom_id], r_s=rs)
                    descriptors.append(temp)

                for lambda_value in Constants.LAMBDA_VALUES:
                    for eta, zeta in product(Constants.ETA_VALUES, Constants.ZETA_VALUES):
                        temp = self._angle_function(zeta=zeta, eta=eta, cosine_theta=cosine_theta,
                                                    distance_ij=distance_ij, distance_jk=distance_jk,
                                                    distance_ik=distance_ik, lambda_value=lambda_value)
                        descriptors.append(temp)
                descriptors_by_atom_id.append(descriptors)
            env_descriptors.append(descriptors_by_atom_id)
        return env_descriptors

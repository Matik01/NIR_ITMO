import concurrent.futures
import argparse
from util.calculate_descriptors import BPS
import os

parser = argparse.ArgumentParser(description='BPS functions calculation tool for PDB files with ligand-protein complex')
parser.add_argument('--pdb_path', type=str, help='Path to PDB files from MD', required=True)
parser.add_argument('--pdb_prefix', type=str, help='Provide MD frames file prefix. Example: "conf" if conf0-1000.pdb ',
                    required=False)



if __name__ == '__main__':
    args = parser.parse_args()
    if args.pdb_prefix is None:
        files = [os.path.join(args.pdb_path, pdb_file) for pdb_file in os.listdir(args.pdb_path) if
                 pdb_file.endswith('.pdb')]
    else:
        files = [os.path.join(args.pdb_path, pdb_file) for pdb_file in os.listdir(args.pdb_path) if
                 pdb_file.endswith('.pdb') and pdb_file.startswith(args.pdb_prefix)]

    for file in files:
        run = BPS(file=file)
        run.main()
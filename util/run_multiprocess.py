from concurrent.futures import ThreadPoolExecutor, as_completed
import h5py
import argparse
from util.calculate_descriptors import BPS
import os

parser = argparse.ArgumentParser(description='BPS functions calculation tool for PDB files with ligand-protein complex')
parser.add_argument('--pdb_path', type=str, help='Path to PDB files from MD', required=True)
parser.add_argument('--pdb_prefix', type=str, help='Provide MD frames file prefix. Example: "conf" if conf0-1000.pdb ',
                    required=False)
parser.add_argument('--n_threads', type=int, help='Number of threads', required=True, default=4)


def save_to_h5(output_file: str, data: list[tuple]):
    with h5py.File(output_file, 'w') as f_out:
        for file_id, array in data:
            f_out.create_dataset(file_id, data=array)


def create_bps_calculations(file_name: str):
    return BPS(file=file_name).main()


if __name__ == '__main__':
    output_file = 'output.h5'
    args = parser.parse_args()
    if args.pdb_prefix is None:
        files = [os.path.join(args.pdb_path, pdb_file) for pdb_file in os.listdir(args.pdb_path) if
                 pdb_file.endswith('.pdb')]
    else:
        files = [os.path.join(args.pdb_path, pdb_file) for pdb_file in os.listdir(args.pdb_path) if
                 pdb_file.endswith('.pdb') and pdb_file.startswith(args.pdb_prefix)]

    with ThreadPoolExecutor(max_workers=args.n_threads) as executor:        # TODO: Better use ProcessPoolExecutor
        futures = {executor.submit(create_bps_calculations, file) for file in files}

        results: list[tuple] = []
        for future in as_completed(futures):
            file_id, result = future.result()
            results.append((file_id, result))

        save_to_h5(output_file, results)

from concurrent.futures import ProcessPoolExecutor, as_completed
import h5py
import argparse
from util.calculate_descriptors import BPS
import os
import traceback
import json

parser = argparse.ArgumentParser(description='BPS functions calculation tool for PDB files with ligand-protein complex')
parser.add_argument('--json_path', type=str, help='Path to PDB files from MD', required=True)


def save_to_h5(output_file: str, data: list[tuple]):
    goal = 'a' if os.path.exists(output_file) else 'w'
    with h5py.File(output_file, goal) as f_out:
        for file_id, array in data:
            f_out.create_dataset(file_id, data=array)


def create_bps_calculations(file_name: str):
    return BPS(file=file_name).main()


def read_json(json_file: str):
    with open(json_file, 'r') as f:
        return json.load(f)


if __name__ == '__main__':
    args = parser.parse_args()
    files = [file.endswith(".json") for temp in os.walk(args.json_path) for file in temp[2]]

    for file in files:
        conf_data = read_json(file)

        if conf_data["pdb_prefix"] is None:
            files = [os.path.join(conf_data["pdb_path"], pdb_file) for pdb_file in os.listdir(conf_data["pdb_path"]) if
                     pdb_file.endswith('.pdb')]
        else:
            files = [os.path.join(conf_data["pdb_path"], pdb_file) for pdb_file in os.listdir(conf_data["pdb_path"]) if
                     pdb_file.endswith('.pdb') and pdb_file.startswith(conf_data["pdb_prefix"])]

        with ProcessPoolExecutor(max_workers=conf_data["n_cpu"]) as executor:
            futures = {executor.submit(create_bps_calculations, file) for file in files}

            results: list[tuple] = []
            for future in as_completed(futures):
                try:
                    file_id, result = future.result()
                    results.append((file_id, result))
                except Exception as e:
                    print(f'Error in file calculation for  {futures[future]}:')
                    traceback.print_exc()

            save_to_h5(conf_data["output_file"], results)

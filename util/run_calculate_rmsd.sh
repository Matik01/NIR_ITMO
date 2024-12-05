#!/bin/bash
#RUN FROM RECEPTOR_GROUP DIRECTORY
for dir in *_package; do
	receptor_id="${dir%_package}"
	echo "Running: $receptor_id"
	python3 ../util/calculate_docking_rmsd.py --workspace_dir "${dir}/workspace/" --ref_ligand "${dir}/data/ligand.pdb" --output_csv "${receptor_id}.csv"
done

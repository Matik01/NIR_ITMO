import os
from rdkit import Chem
from rdkit.Chem import rdMolAlign
import pandas as pd
import re
import glob
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import argparse
import sys


parser = argparse.ArgumentParser(
    prog='calcDockRmsd'
)

parser.add_argument('--workspace_dir')
parser.add_argument('--ref_ligand')
parser.add_argument('--output_csv')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

rootdir = args.workspace_dir
directories = []
for dir, subdir, file in os.walk(rootdir):
    if len(dir.split(sep='/')) == 4:
        directories.append(dir)


id_best_rmsd = dict()
ref_ligand = glob.glob(args.ref_ligand)
ref_mol = Chem.MolFromPDBFile(ref_ligand[0])
ref_smarts = Chem.MolToSmarts(ref_mol)

for dir in directories:
    exp_id = dir.split(sep='/')[2].split(sep="@")[0]
    exp_directory = dir.split(sep='/')[2]
    id_best_rmsd[exp_id] = [0, 9999] #pose/rms values
    pdb_files = []
    main_dir = None

    for dir_exp, subdir, files in os.walk(dir):
        if main_dir is None:
            main_dir = dir_exp
        pdb_files.extend([file for file in files if "pml.pdb" in file])
    
    for conf_pdb in pdb_files:
        mols = []
        mol = Chem.MolFromPDBFile(f'{dir}/{conf_pdb}', removeHs=True)
        mols.append(ref_mol)
        mols.append(mol)
        mcs = rdFMCS.FindMCS(mols)
        pattern = Chem.MolFromSmarts(mcs.smartsString)
        ref = mols[0]
        refMatch = ref.GetSubstructMatch(pattern)
        molMatch = mol.GetSubstructMatch(pattern)
        rms = Chem.rdMolAlign.CalcRMS(mol, ref_mol, map=[list(zip(molMatch, refMatch))])
        if rms < id_best_rmsd.get(exp_id)[1]:
            if id_best_rmsd.get(exp_id)[1] == 9999:
                pass
            else:
                os.remove(f'{main_dir}/{exp_directory}.{id_best_rmsd.get(exp_id)[0]}.pml.pdb')
            id_best_rmsd[exp_id] = [conf_pdb.split(sep='.')[1], rms]
        elif rms > id_best_rmsd.get(exp_id)[1]: 
            os.remove(main_dir + "/" + conf_pdb)

df = pd.DataFrame(
    [(key, value[0], value[1]) for key, value in id_best_rmsd.items()], 
    columns=['id', 'pose', 'rmsd_value']
    )
df.to_csv(args.output_csv)

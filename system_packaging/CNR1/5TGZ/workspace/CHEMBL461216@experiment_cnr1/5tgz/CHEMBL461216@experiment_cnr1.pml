load CHEMBL461216@experiment_cnr1.dok.pdbqt
alter resn UNL, resn="SUS"
python
for i in range(1,cmd.count_states()+1): cmd.save(f'CHEMBL461216@experiment_cnr1.{i}.pml.pdb', selection="resn SUS", state=i)
python end
quit

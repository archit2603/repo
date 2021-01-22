import pandas as pd
import sys
import os

dataset = sys.argv[1]

# load smiles
compound_iso_smiles = []
for smi in os.listdir('data/'+dataset+'/smiles'):
    
    with open('data/'+dataset+'/smiles/'+smi, 'r') as file:
        data = file.read()

    for comp in data.split('\n'):
        if(comp.split(' ')[0]!='' and comp.split(' ')[0]!='SMILES' ):
            compound_iso_smiles.append(comp.split(' ')[0])
compound_iso_smiles = list(set(compound_iso_smiles))


# laod targets
pdbs = []
pdbs_seqs = []
for seq in os.listdir('data/'+dataset+'/sequences'):
    
    with open('data/'+dataset+'/sequences/'+seq, 'r') as file:
        data = file.read()

    pdbs.append(data.split('\n')[0])
    pdbs_seqs.append(data.split('\n')[1])

# create dataset
compound_iso_smiles_temp = []
pdbs_temp = []
pdbs_seqs_temp = []
num_targets = len(pdbs)
num_smiles = len(compound_iso_smiles)
temp = []
for i in range(num_targets):
    compound_iso_smiles_temp += compound_iso_smiles
    
    temp = [pdbs[i][1:]]
    temp = temp * num_smiles
    pdbs_temp += temp
    
    temp = [pdbs_seqs[i]]
    temp = temp * num_smiles
    pdbs_seqs_temp += temp


pdbs = pdbs_temp
pdbs_seqs = pdbs_seqs_temp
compound_iso_smiles = compound_iso_smiles_temp
del(temp)
del(compound_iso_smiles_temp)
del(pdbs_temp)
del(pdbs_seqs_temp)


dict = {'compound_iso_smiles': compound_iso_smiles, 'target_name': pdbs, 'target_sequence': pdbs_seqs}
df = pd.DataFrame(dict)
df.to_csv('data/'+dataset+'/split/'+dataset+'.csv')

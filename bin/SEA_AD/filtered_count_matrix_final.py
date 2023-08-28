##To randomly sample cells of an AD stage

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.sparse import csr_matrix
import gzip
import itertools
import random
import sys
import math
import json


mtg_rna = h5py.File('../../AD_data/SEA_AD/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2023-05-05.h5ad')

gene_id = mtg_rna.get('var').get('gene_ids')
gene_id = pd.DataFrame(data=gene_id,columns=['gids'])
gene_id=gene_id.iloc[:,0].apply(lambda s: s.decode('utf-8'))

X_data = mtg_rna.get('X').get('data')
X_indices = mtg_rna.get('X').get('indices')
X_indptr = mtg_rna.get('X').get('indptr')

braak = mtg_rna.get('obs').get('Braak')
braak = pd.DataFrame(data=braak)
ref = list(np.where(braak.iloc[:,0]==0)[0])
br_6 = list(np.where(braak.iloc[:,0]==6)[0])

num_samples = 50000
ref_samp = sorted(random.sample(ref, num_samples))
br6_samp = sorted(random.sample(br_6, num_samples))

global_indices_ref = {}
for i in ref_samp:
    global_indices_ref[int(i)] = list(range(int(X_indptr[i]),int(X_indptr[i+1])))

global_indices_br6 = {}
for i in br6_samp:
    global_indices_br6[int(i)] = list(range(int(X_indptr[i]),int(X_indptr[i+1])))
global_indices_ref.update(global_indices_br6)

with open('../../results/SEA_AD/global_indices_final.json', 'w') as convert_file:
     convert_file.write(json.dumps(global_indices_ref))

meta = len(global_indices_br6)*['ref'] + len(global_indices_br6)*['br6']
with open('../../results/SEA_AD/global_inidces_meta_final.txt', 'w') as f:
    for line in meta:
        f.write(f"{line}\n")

global_indices_final = global_indices_ref

#for count,col in enumerate(gene_id):
    #if count==len(gene_id)-1:
        #print(col)                                                                                                                                                                                        
    #else:                                                                                                                                                                                                 
        #print(col,end='\t')

counts = pd.DataFrame(columns=[list(range(len(gene_id)))])
for cell in global_indices_final.keys():
    counts.loc[cell,X_indices[global_indices_final[cell]]] = X_data[global_indices_final[cell]]


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 14:36:06 2023

merge results from dompair-eval-main into one big dict.

@author: Shaun Kandathil, 2023
"""
import os, sys
import pickle

import numpy as np
from tqdm import tqdm

iv_data_dir = "/cluster/project9/afdb_domain_ext/results/dompair-eval-outputs-v2"
out_dir = iv_data_dir

out_file = os.path.join(out_dir, "all_interaction_vectors_raw_v2.pkl")

master = dict()
os.chdir(iv_data_dir)
chunk_count=0

for i in tqdm(range(1, 5133)):
    fname = 'dompair_interaction_vectors_chunk' + str(i).zfill(6) + '_with_choppings.pkl'
    try:
        p = pickle.load(open(fname, 'rb'))
        for k in p.keys():
            if k not in master.keys():
                master[k] = p[k] # also sets ref_domain_pair (the ref is a single domain, not a pair; correct info is in refdom_d anyway)
                del master[k]['ref_domain_pair']
            else:
                master[k]['vectors'] = np.concatenate([ master[k]['vectors'], p[k]['vectors'] ])
                #master[k]['aligned_domain_pairs'].extend(p[k]['aligned_domain_pairs'])
                for y in ['aligned_domain_pairs', 'choppings', 'pae_score']:
                    master[k][y].extend(p[k][y])
        chunk_count+=1
    except FileNotFoundError:
        print("Could not find", fname, file=sys.stderr)
        
    
print("Processed", chunk_count, "chunks, npify pae_score...")

for k in master.keys():
	master[k]['pae_score'] = np.array(master[k]["pae_score"])

pickle.dump(master, file=open(out_file, 'wb'))
print("All vectors written to", out_file)

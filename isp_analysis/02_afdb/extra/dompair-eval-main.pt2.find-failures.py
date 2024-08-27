#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# iterate thru all output chunks and mark as failure if:
# 1. the file doesn't exist, or is missing subkey 'choppings' or its length is off (this requires update_iv_chunks_with_choppings.py to be rerun)
# 2. the file does not have subkey 'pae_score' (rerun either variant of dompair-eval-main.pt2.py)
# 3. length of the list under 'pae_score' != length of list under 'aligned_domain_pairs', or subkey 'pae_score' missing (remedy as in 2)

###############################################################################
## IMPORTS

import sys,os
import pickle

from tqdm import tqdm

###############################################################################
## GLOBALS

basedir = os.path.join(os.environ['HOME'], 'afdb_domain', 'tools', 'dompair-eval-afdb')

out_dir = "/cluster/project9/afdb_domain_ext/results/dompair-eval-outputs-v2"
os.makedirs(out_dir, exist_ok=True)

# separator used to join dbids and domnames
sep = '-'
domname_sep = ':'

# string in pdb filenames
pdb_mid_str = '-model_v4_'

# input/output filename base
iv_out_f_base = '/'.join([out_dir,'dompair_interaction_vectors_chunk'])

failure_outf = "failures.pt2"
failure_outf2 = "failures.pt2.noinput"

count1 = 0
count2 = 0

total = 5132

with open(failure_outf, 'w') as f, open(failure_outf2, 'w') as f2:
    for chunk_id_ in tqdm(range(1,total+1)):
    
        chunk_id = str(chunk_id_).zfill(6)
        iv_out_f2 = iv_out_f_base + chunk_id + "_with_choppings.pkl"
        if os.path.exists(iv_out_f2):
            iv = pickle.load(open(iv_out_f2, 'rb'))
        else:
            print(chunk_id_, file=f2)
            count2 +=1
            continue
        
        for isp in iv.keys():
            if 'pae_score' not in iv[isp].keys() or len(iv[isp]['pae_score']) != len(iv[isp]['aligned_domain_pairs']):
                print(chunk_id_, file=f)
                count1+=1
                break


if count1 > 0:
    print("Class 1 failures in",failure_outf, "(rerun dompair-eval-main.pt2.py):", count1)
if count2 > 0:
    print("Class 2 failures in",failure_outf2,"(rerun update_iv_chunks_with_choppings.py):", count2)
if count1 == count2 == 0:
    print('No failures!')
    os.remove(failure_outf)
    os.remove(failure_outf2)

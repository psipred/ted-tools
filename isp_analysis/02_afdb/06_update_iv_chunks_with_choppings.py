#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Add choppings to dompair eval output pickles, we need these to process PAE data and filter.

@author: Shaun M Kandathil

"""
###############################################################################
## IMPORTS

import sys,os
#import shutil
import pickle
#import subprocess
#from itertools import combinations, chain
#from copy import deepcopy

import numpy as np
from tqdm import tqdm
# from Bio.PDB import PDBParser, NeighborSearch
# from Bio.PDB.Structure import Structure

#import tmalign
import pae_util

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


###############################################################################
## FUNCTIONS


###############################################################################
## MAIN

big_index_f = "afdb_ids_with_isps.lst.lookupbyid.pkl"

print("Loading big index with choppings...")
idx = pickle.load(open(big_index_f, 'rb'))
print("Done.")

print("Commence update...")

total = 5132

for chunk_id_ in tqdm(range(1,total+1)):

    chunk_id = str(chunk_id_).zfill(6)

    iv_out_f = iv_out_f_base + chunk_id + ".pkl"
    iv_out_f2 = iv_out_f_base + chunk_id + "_with_choppings.pkl"

    iv = pickle.load(open(iv_out_f, 'rb'))

    for isp in iv.keys():
        iv[isp]["choppings"] = list()

        for adp in iv[isp]["aligned_domain_pairs"]:

            s = adp.split(domname_sep)
            afdb_id = adp.split(pdb_mid_str)[0]

            # get the choppings
            choppings = [ idx[afdb_id]['choppings'][x] for x in s ]
            iv[isp]["choppings"].append(domname_sep.join(choppings))

    assert len(iv[isp]["choppings"]) == len(iv[isp]["aligned_domain_pairs"])
    pickle.dump(iv, file=open(iv_out_f2, 'wb'))


print("All done.")

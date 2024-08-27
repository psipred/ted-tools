#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 11:02:38 2023
# Add pae_score to dompair pkls.
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
chunk_id = sys.argv[1] #"000001" # or just '1' #
chunk_id = chunk_id.zfill(6)

iv_out_f = iv_out_f_base + chunk_id + "_with_choppings.pkl"

#big_index_f = "afdb_ids_with_isps.lst.lookupbyid.pkl"

print("Loading result chunk", chunk_id)
iv = pickle.load(open(iv_out_f, 'rb'))
print("Done.")

# print("Loading big index with choppings...")
# idx = pickle.load(open(big_index_f, 'rb'))
# print("Done.")

print("Commence filtering...")

for isp in iv.keys():
    iv[isp]["pae_score"] = list()

    for i in range(len(iv[isp]["aligned_domain_pairs"])):
        adp = iv[isp]["aligned_domain_pairs"][i]
        s = adp.split(domname_sep)
        afdb_id = adp.split(pdb_mid_str)[0]

        # get the choppings
        # chopping_d1, chopping_d2 = [ idx[afdb_id]['choppings'][x] for x in s ]
        chopping_d1, chopping_d2 = iv[isp]["choppings"][i].split(domname_sep)
        # convert to resranges
        res_d1 = pae_util.chopping_to_res(chopping_d1)
        res_d2 = pae_util.chopping_to_res(chopping_d2)
        
        # get the pae matrix
        try:
            pae = pae_util.get_pae_from_afdb_id(afdb_id)

            # get submatrix and median
            mdn_pae = pae_util.interdom_pae_summary(pae, res_d1, res_d2, both=True)
            iv[isp]["pae_score"].append(mdn_pae)
        except Exception as e:
            # we don't write the output pkl until the loops complete,
            # so this is redundant, but we may want to do sth different in the future
            print(e)
            sys.exit(1)

assert len(iv[isp]["pae_score"]) == len(iv[isp]["aligned_domain_pairs"])

print("All ops complete; overwriting", iv_out_f)
pickle.dump(iv, file=open(iv_out_f, 'wb'))

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
add PAE score to pkls from local tar files
@author: Shaun M Kandathil
"""
###############################################################################
## IMPORTS

import sys,os
import shutil
import pickle
import subprocess
import json
import gzip
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

#chaindomdic (lookup-dict) prefix and suffix. use these to get proteomes
cdd_prefix = os.path.join("/cluster/project9/afdb_domain_ext/results/dompair-eval-inputs-chunks",'lookup-dict-chunk')
cdd_suffix = '.pkl'

out_dir = "/cluster/project9/afdb_domain_ext/results/dompair-eval-outputs-v2"
os.makedirs(out_dir, exist_ok=True)

# I say tgz but actually is tar
tgz_dir = "/cluster/ref3/alphafolddb/proteomes"
tgz_prefix = "proteome-tax_id-"
tgz_suffix = "_v4.tar"
# proteome-tax_id-53326-6_v4.tar

#  AF-A0A023HHL0-F1-predicted_aligned_error_v4.json.gz
pae_file_suffix = "-predicted_aligned_error_v4.json.gz"

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
chunk_id = sys.argv[1] #"000001" # or just '1'
chunk_id = chunk_id.zfill(6)

tmpdir = sys.argv[2]
os.makedirs(tmpdir, exist_ok=True)

cdd_fname = cdd_prefix + chunk_id + cdd_suffix

print('processing', cdd_fname)
print('Load input chunk and invert...')

chaindomdic = pickle.load(open(cdd_fname,'rb'))
cdd_inv = dict()
for proteome in chaindomdic.keys():
    for chain in chaindomdic[proteome].keys():
        cdd_inv[chain] = proteome

print('chaindomdic inversion complete.')


iv_out_f = iv_out_f_base + chunk_id + "_with_choppings.pkl"

print("Loading result chunk", chunk_id)
iv = pickle.load(open(iv_out_f, 'rb'))
print("Done.")

print("Commence filtering...")

for isp in iv.keys():
    iv[isp]["pae_score"] = list()

    for i in range(len(iv[isp]["aligned_domain_pairs"])):
        adp = iv[isp]["aligned_domain_pairs"][i]
        s = adp.split(domname_sep)
        afdb_id = adp.split(pdb_mid_str)[0]

        proteome = cdd_inv[afdb_id]

        # get tar file and extract pae json.gz to tmpdir if not already there
        tgz_base     = tgz_prefix + proteome + tgz_suffix
        tgz_path     = tgz_dir + '/' + tgz_base
        tgz_f_in_tmp = tmpdir  + "/" + tgz_base

        if not os.path.exists(tgz_f_in_tmp):
            shutil.copy(tgz_path, tgz_f_in_tmp)

        pae_f = afdb_id + pae_file_suffix
        
        if not os.path.exists(pae_f):
            tarcmd = ['tar', '-C', tmpdir, '-xf', tgz_f_in_tmp, pae_f]
            tarproc = subprocess.run(args=tarcmd, capture_output=False)
            if tarproc.returncode != 0:
                print('untar error:', file=sys.stderr)
                # print('in dir '+tmpdir, file=sys.stderr)
                print('command:',tarcmd, file=sys.stderr)
                sys.exit(1)

        # get the choppings
        # chopping_d1, chopping_d2 = [ idx[afdb_id]['choppings'][x] for x in s ]
        chopping_d1, chopping_d2 = iv[isp]["choppings"][i].split(domname_sep)

        # convert to resranges
        res_d1 = pae_util.chopping_to_res(chopping_d1)
        res_d2 = pae_util.chopping_to_res(chopping_d2)
        
        # get the pae matrix
        # pae = pae_util.get_pae_from_afdb_id(afdb_id)
        pae_dict = json.load(gzip.open(os.path.join(tmpdir, pae_f), 'rt'))[0]
        pae = np.asarray(pae_dict['predicted_aligned_error'])

        # get submatrix and median
        mdn_pae = pae_util.interdom_pae_summary(pae, res_d1, res_d2, both=True)
        iv[isp]["pae_score"].append(mdn_pae)

assert len(iv[isp]["pae_score"]) == len(iv[isp]["aligned_domain_pairs"])

print("All ops complete; overwriting", iv_out_f)
pickle.dump(iv, file=open(iv_out_f, 'wb'))

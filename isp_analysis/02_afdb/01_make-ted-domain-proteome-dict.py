#!/usr/bin/env python3

# make dicts with proteome files as keys, AFDB ids and TED domain ids as values
# Shaun Kandathil, 2023
##########################################################################################
#### IMPORTS

import os
import pickle
import time
from copy import deepcopy

import numpy as np

##########################################################################################
#### GLOBALS

#os.chdir("/cluster/project9/afdb_domain_ext/results")

# string that joins AFDB entry name and TED domain name (e.g. 'AF-A0A016WTH9-F1-model_v4_TED01')
joiner = '-model_v4_'

#string prefix and suffix for proteome filenames
pprefix="proteome-tax_id-"
psuffix="_v4.consensus_domains"

#source file with ONLY multi-domain entries
# columns in srcfile:
# domain_id domain_md5 consensus_level chopping domain_nres n_segments average_plddt n_ss n_helix n_strand n_helix+n_strand n_turn proteome-domain_filename ted_cath_id
srcfile = "h_combined_hits.domain_summary.sort.cath.multidomain.tsv"

# SMK srcfile must be in the same format as above, but may contain a subset of the lines in it.
# For now, process the entire file as the `proteome-domain-filename` field is not sorted.
# Then once all lines have been read, we spit out multiple files, each containing `chunksize` chains.

pklfile_base = "dompair-eval-inputs-chunks/lookup-dict-chunk"

os.makedirs('dompair-eval-inputs', exist_ok=True)
os.makedirs('dompair-eval-inputs-chunks', exist_ok=True)

# check that pklfile_base is writable...
with open(pklfile_base, 'w') as f:
    f.write('foo')
f.close()
os.remove(pklfile_base) # because every inode counts!

consensus_level_int = {'high':2, 'medium':1, 'low':0}
high_consens = consensus_level_int['high']

# number of *chains* in each output pkl.
# only a small fraction of domain pairs are in contact (4% in a sample of 1000 chains)
chunksize = 10000

##########################################################################################
#### FUNCTIONS

def npify(elem):
    for k in elem.keys():
        for k2 in elem[k].keys():
            elem[k][k2] = np.array(elem[k][k2])
    return elem


def mutate_dict(f,d):
    for k, v in d.items():
        d[k] = f(v)

##########################################################################################
#### MAIN
print('Commence reading...')

start_time = time.time()
global_start_time = start_time
d = dict()
refdom_d = dict()
best_consens = dict()
counter=0


with open(srcfile, 'r') as f:
    for line in f:
        data = line.rstrip('\n').split()
        # assert len(data) == 14
        pfilename = data[-2] # e.g. proteome-tax_id-1445577-1_v4.consensus_domains
        shortp = pfilename.removeprefix(pprefix).removesuffix(psuffix)
        ted_id = data[-1] # e.g. 3.40.50.720
        full_id = data[0].split(joiner)
        plddt = data[6]

        chain_id = full_id[0] # e.g. 'AF-A0A016WTH9-F1'
        domid = full_id[1]  # e.g. 'TED01'

        consensus_level = consensus_level_int[data[2]]

        if shortp not in d.keys():
            d[shortp] = dict()
            # This condition is implied by the above condition being true, so omit:
            # if chain_id not in d[shortp].keys():
            d[shortp][chain_id] = dict()
            d[shortp][chain_id]['domids'] = [domid]
#            d[shortp][chain_id]['consensus'] = [consensus_level]
            d[shortp][chain_id]['tedids'] = [ted_id]
            d[shortp][chain_id]['plddts'] = [plddt]
            counter +=1
        else:
            if chain_id not in d[shortp].keys():
                d[shortp][chain_id] = dict()
                d[shortp][chain_id]['domids'] = [domid]
#                d[shortp][chain_id]['consensus'] = [consensus_level]
                d[shortp][chain_id]['tedids'] = [ted_id]
                d[shortp][chain_id]['plddts'] = [plddt]
            else:
                d[shortp][chain_id]['domids'].append(domid)
#                d[shortp][chain_id]['consensus'].append(consensus_level)
                d[shortp][chain_id]['tedids'].append(ted_id)
                d[shortp][chain_id]['plddts'].append(plddt)

        if ted_id not in refdom_d.keys():
            refdom_d[ted_id] = dict()
            refdom_d[ted_id]['taxid'] = shortp
            refdom_d[ted_id]['chain'] = data[0]
            refdom_d[ted_id]['plddt'] = plddt
            # refdom_d[ted_id]['cons_lvl'] = consensus_level
        else:
            if plddt > refdom_d[ted_id]['plddt']:
                refdom_d[ted_id]['taxid'] = shortp
                refdom_d[ted_id]['chain'] = data[0]
                refdom_d[ted_id]['plddt'] = plddt
                # refdom_d[ted_id]['cons_lvl'] = consensus_level

        if ted_id not in best_consens.keys():
            best_consens[ted_id] = consensus_level
        else:
            if consensus_level > best_consens[ted_id]:
                best_consens[ted_id] = consensus_level


elapsed = time.time() - start_time

print("Processed {} keys in {} sec, npify...".format(counter, elapsed))

start_time = time.time()

mutate_dict(npify, d)

elapsed = time.time() - start_time

print('Npified in {} sec, create sub-index files...'.format(elapsed))

start_time = time.time()


chunk=1
d_out = dict()
chaincount=0
for k in d.keys():
    for chain in d[k].keys():
        # k is a proteome id, chain is a chainid
        if k not in d_out.keys():
            d_out[k] = dict()

        # all domains for any chain are always in the same output dict.
        d_out[k][chain] = deepcopy(d[k][chain])
        chaincount+=1
        if chaincount >= chunksize:
            # open new chunk file
            chunk_f = pklfile_base+str(chunk).zfill(6)+".pkl"
            pickle.dump(d_out, open(chunk_f, 'wb'))
            chunk+=1
            chaincount=0
            d_out=dict()

if len(d_out) != 0:
    # open new chunk file
    chunk_f = pklfile_base+str(chunk).zfill(6)+".pkl"
    pickle.dump(d_out, open(chunk_f, 'wb'))

# write out refdom_d and best_consens

pickle.dump(refdom_d, open('dompair-eval-inputs/ref_domains_by_tedid.pkl','wb'))
pickle.dump(best_consens, open('dompair-eval-inputs/best_consens.pkl', 'wb'))

elapsed = time.time() - start_time

print('Last chunk number:', chunk)

print("Files written in {} sec.".format(elapsed))

total_time = time.time() - global_start_time

print("All done; total time {} sec.".format(total_time))

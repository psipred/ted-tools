#!/usr/bin/env python3

# create a data structure that lets us look up which cath IDs are in which pdb chain.

# import sys, os
import pickle
import numpy as np

fn = 'cath-domain-list.txt'

def npify(elem):
    for k in elem.keys():
        elem[k] = np.array(elem[k])
    return elem

master={}

# some terminology:
# domid: e.g. 1acfA00
# dbid/cathid: e.g. 1.10.5.11
# chain: e.g. 1acfA

with open(fn, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue

        arr = line.strip().split()

        # if arr[0][5:7] == '00':
        #     continue # single-domain chain
        domid = arr[0]

        chain = arr[0][0:5] # CATH4.3 doesn't seem to have chain ids with more than 1 char...

        cathid = '.'.join(arr[1:5])

        if chain in master.keys():
            master[chain]['domids'].append(domid)
            master[chain]['cathids'].append(cathid)
        else:
            master[chain]={}
            master[chain]['domids'] = [domid]
            master[chain]['cathids'] = [cathid]


master = {k:npify(v) for k,v in master.items()}

pickle.dump(master, open(fn+'.withdomids.pkl', 'wb'))

# also make a master list of domains. Any/all domains we will need to process will be a strict subset of these.

l = []

print("Reading doms...")

for _, v in master.items():
    l.extend(list(v['cathids']))
print('creating set of unique doms...')

domlist = set(l)

print(f'There are {len(domlist)} unique CATH superfamilies.')

pickle.dump(domlist, open(fn+'.uniquecathids-set.pkl', 'wb'))

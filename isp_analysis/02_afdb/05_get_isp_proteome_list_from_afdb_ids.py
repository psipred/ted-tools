#!/usr/bin/env python3
import os
import pickle
# from copy import deepcopy

inputf = "afdb_ids_with_isps.lst"
outputf = inputf + '.pkl'

all_ids = None

# create set of required afdb ids for fast lookup
if os.path.exists(outputf):
    all_ids = pickle.load(open(outputf, 'rb'))
    print('Loaded ids from'+outputf)
else:
    all_ids = list()
    print("Reading IDs from",inputf,"and storing as set...")
    with open(inputf, 'r') as f:
        for line in f:
            all_ids.append(line.rstrip())

    all_ids = set(all_ids)
    
    pickle.dump(all_ids, file=open(outputf, 'wb'))

print("gather...")
# search against zip index and print(?)
resultbyproteome = dict()
resultbyid = dict()

# this file has data for all 214M AFDB targets in three columns: AFDB_id, seq_MD5, proteome_id
with open(os.environ['HOME'] + '/afdb_domain/datasets/index/afdb_214m_targets.zip_index', 'r') as f:
    for line in f:
        s = line.strip().split()
        afdb_id = s[0].split('-model_v4')[0]
        proteome = s[2]
        if afdb_id in all_ids:
            if proteome not in resultbyproteome.keys():
                resultbyproteome[proteome] = dict()
            resultbyproteome[proteome][afdb_id] = dict() # we add info about domains below

            assert afdb_id not in resultbyid.keys()

            resultbyid[afdb_id] = {"proteome":proteome, "choppings": dict()}


# also get choppings
with open("h_combined_hits.domain_summary.sort.cath.multidomain.tsv", 'r') as f:
    for line in f:
        s = line.rstrip().split()
        dom = s[0]
        chopping = s[3]
        proteome = s[12].removesuffix(".consensus_domains")
        afdb_id = dom.split("-model_v4")[0]

        if proteome in resultbyproteome.keys():
            if afdb_id in resultbyproteome[proteome].keys():
                resultbyproteome[proteome][afdb_id][dom] = chopping
                resultbyid[afdb_id]["choppings"][dom] = chopping


print("write index by id")
pickle.dump(resultbyid, file=open(inputf+".lookupbyid.pkl", 'wb'))

del resultbyid

print("write index by proteome")
pickle.dump(resultbyproteome, file=open(inputf+".lookupbyproteome.pkl", 'wb'))

lstf = "isp_proteomes.lst"
print("Write", lstf)
with open(lstf, 'w') as f:
    for k in resultbyproteome.keys():
        print(k, file=f)
f.close()

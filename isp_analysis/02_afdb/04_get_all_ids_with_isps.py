#!/usr/bin/env python3

import pickle

def get_afdbs(adp):
    return set([x.split('-model_v4')[0] for x in adp['aligned_domain_pairs']])


d = pickle.load(open('all_interaction_vectors_raw_v2.pkl','rb'))
outf = 'afdb_ids_with_isps.lst'

l = [get_afdbs(x) for x in d.values()]

all_afdb_ids = set().union(*l)

with open(outf, 'w') as f:
    for i in all_afdb_ids:
        print(i, file=f)

print(len(all_afdb_ids),' IDs written to', outf)

pickle.dump(all_afdb_ids, open(outf+".pkl", 'wb'))
print("Also created", outf+".pkl")

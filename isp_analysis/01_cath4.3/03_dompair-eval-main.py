#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 11:02:38 2023

@author: Shaun M Kandathil
"""
###############################################################################
## IMPORTS

import sys,os
import shutil
import pickle
from itertools import combinations, chain
from copy import deepcopy

import numpy as np
from Bio.PDB import PDBParser, NeighborSearch
from tqdm import tqdm

import tmalign

###############################################################################
## GLOBALS

db = 'cath'

basedir = os.path.join(os.environ['HOME'], 'afdb_domain/tools/dompair-eval-cath4.3')

# output filename
iv_out_f = os.path.join(basedir,'dompair_interaction_vector_data.pkl')

# Path containing domain PDB files (filenames do not have '.pdb' extension, e.g. 1acfA00
dom_pdb_dir = os.path.join(os.environ['HOME'], 'cath-dompdb')

tmpdir = '/dev/shm'

# atom distance thresholds in Angstrom for Cbeta and any-heavy-atom contacts
contact_params = {'cbeta': 8.0, 'heavy': 3.0}

# contact definition based on C-beta atom ('cbeta') or any heavy atom ('heavy')
cont_type = 'cbeta'

# min number of contacts between 2 domains to say that they are in contact
min_contacts = 3

# separator used to join dbids and domnames
sep = '-'
domname_sep = ':'

# path to tmalign executable
tmalign_bin = os.path.join(os.environ['HOME'], "software_staging", "tmalign", "tmalign")

# copy executable to faster, more reliable storage
username = os.environ['LOGNAME']
if os.path.exists("/scratch0"):
    os.makedirs("/scratch0/"+username, exist_ok=True)
    shutil.copy2(tmalign_bin, "/scratch0/"+username+"/tmalign")
    tmalign_bin = "/scratch0/"+username+"/tmalign"

# tmalign args (NEVER set -m or -outfmt)
tmalign_args = "-ter 0"

# Enable warnings
warn = True

# if `warn`, warn when TM-align finds a TM-score of this much or less when
#   aligning two domains from the same superfamily
tmscore_warn_thresh = 0.4

###############################################################################
## FUNCTIONS

def dompair_to_id(doms, sep=sep):
    doms=list(doms)
    doms.sort()
    return sep.join(doms)


def id_to_dompair(id_, sep=sep):
    s = id_.split(sep)
    if len(s) != 2:
        raise AssertionError (id_+" did not split into 2 using '"+sep+"'." )
        sys.exit(1)
    return s


def sort_list(t):
    # take an iterable, convert to list and then sort it
    l = list(t)
    l.sort()
    return l


# mutate a dict in-place to save memory (hopefully)
def mutate_dict(f,d):
    for k, v in d.items():
        d[k] = f(v)


# def make_combinations(d:pd.DataFrame):
#     nrow = len(d.index)
#     if nrow == 2:
#         return [d]

#     combns = combinations(range(nrow), 2) # is an iterator

#     return [d.loc[c,:] for c in combns] # each c is a 2-tuple of row indices


def get_atoms_of_interest(chain1, chain2, mode='cbeta'):
    aoi = []
    if mode == "cbeta":
        for r in chain(chain1.get_residues(), chain2.get_residues()):
            a = 'CA' if r.get_resname() == 'GLY' else 'CB'
            try:
                aoi.append(r[a])
            except KeyError: # fall back to CA if CB atom is missing
                try:
                    aoi.append(r['CA'])
                except KeyError:
                    pass  # if that fails as well then just skip this residue

    elif mode == "heavy":
        for r in chain(chain1.get_residues(), chain2.get_residues()):
            for a in r:
                if a.element == 'H':
                    continue
                aoi.append(r[a])
    else:
        print('In get_atoms_of_interest(): "mode" must be either "cbeta" or "heavy"', file=sys.stderr)
        sys.exit(1)

    return aoi


def domains_in_contact(dom1, dom2, min_contacts=3, cont_type='cbeta'):
    global contact_params
    aoi = get_atoms_of_interest(dom1, dom2, cont_type)
    ns = NeighborSearch(aoi)
    nbrs = ns.search_all(contact_params[cont_type])
    obs_cont = 0

    for nbr in nbrs:
        if nbr[0].full_id[0] != nbr[1].full_id[0]:
            obs_cont += 1
        if obs_cont >= min_contacts:
            return True

    return False
    # return [b for b in map(lambda a: (a[0].full_id[3][1], a[1].full_id[3][1]), nbrlist)]


def txtfile_to_list(fname):
    """
    Take lines from a text file and place in a list
    Could do with readlines() but that could cause problems on large files
    """
    l=[]
    with open(fname, 'r') as f:
        for line in f:
            l.append(line.strip())
        f.close()
    return l

###############################################################################
## MAIN

if db == 'cath':
    domset = pickle.load(open(os.path.join(basedir, "cath-domain-list.txt.uniquecathids-set.pkl"),'rb'))
    chaindomdic = pickle.load(open(os.path.join(basedir, "cath-domain-list.txt.withdomids.pkl"),'rb'))

    dompair_info_f = "cath-dompair-info.pkl"

    dbids = 'cathids'
elif db == 'ecod':
    raise NotImplementedError("You'll need to write more code to deal with ECOD instead of CATH.")
#    domset = pickle.load(open(os.path.join(basedir, "ecod.latest.domains.txt.unique-ecod-ids-set.pkl"),'rb'))
#    chaindomdic = pickle.load(open(os.path.join(basedir, "ecod.latest.domains.txt.withdomids.pkl"),'rb'))
#    dbids = 'ecodids'
else:
    sys.exit('invalid db - only CATH is currently supported.')



print("Total number of chains:", len(chaindomdic.keys()))

# first filter the chaindomdic to only retain chains with 2 or more domains.
chaindomdic = {k:v for k,v in chaindomdic.items() if len(v["domids"]) > 1}

chainlist = chaindomdic.keys() ## make cmdline param to work on a subset


print("Total number of multi-domain chains:", len(chaindomdic.keys()))

# this is a dict of simply which cathid pairs are *present* in which chains.
# They are not filtered based on being in contact.
pairinfo = dict()

if os.path.exists(dompair_info_f):
    pairinfo = pickle.load(open(dompair_info_f, 'rb'))
    print('Loaded pairinfo pkl with', len(pairinfo.keys()), "unique pairs of", dbids)
else:
    print('pairinfo pkl not found; create from scratch...')
    # possible pairinfo structure:
    # lvl1 key: lexically sorted cathid pair as string (e.g. "1.1.1.1-2.2.2.2")
    # lvl2 key: "chainlist": list() or set(); chains with this pair in it
    # lvl2 key: "count": int; number of chains in pairinfo[cathid_pair]["chainlist"]
    for k in chaindomdic.keys():

        dbid_list = chaindomdic[k][dbids]
        combns = list(map(sort_list, combinations(dbid_list, 2)))
        dompair_ids = set(list(map(dompair_to_id, combns)))
        for i in dompair_ids:
            if not i in pairinfo.keys():
                pairinfo[i] = dict()
                pairinfo[i]['chainlist'] = set()
            pairinfo[i]['chainlist'].add(k)


        #     for domid in chaindomdic[k]['domids']:
        #         dom_pdb_f = os.path.join(dom_pdb_dir, domid)
        #         if not os.path.exists(dom_pdb_f) or os.path.getsize(dom_pdb_f) == 0:
        #             print(k)
        #         i=i+1
        #         if i%1000 == 0:
        #             print("completed", i)

    pickle.dump(pairinfo, open(dompair_info_f, 'wb'))
    print('pairinfo pkl dumped to', dompair_info_f)


result = dict()
chaincount=1
nchain = len(chainlist)

for k in tqdm(chainlist, total=nchain):
    #print(k, "(", chaincount, "/", nchain, ")")
    chaincount+=1

    parser = PDBParser(QUIET=True)
    chain_dom_structs = dict() # name-addressable list

    # load (all) domain pdbs for the chain - I/O is slow
    domnames  = chaindomdic[k]['domids']
    dom_dbids = chaindomdic[k][dbids]

    dom_pdb_fnames = [os.path.join(dom_pdb_dir, a) for a in domnames]

    for f,d in zip(dom_pdb_fnames, domnames):
        f_in_tmp = os.path.join(tmpdir, d)
        if not os.path.exists(f_in_tmp):
            shutil.copy(f, f_in_tmp)
        # read
        chain_dom_structs[d] = parser.get_structure(d, f_in_tmp)

    ndom = len(domnames)

    # find contacting pairs of domains
    combns = combinations(list(range(ndom)), 2)

    for pair in combns:
        i,j = pair # indices into each of the 2 arrays in chaindomdic[k]

        pair_domnames  =  domnames[pair,]  # [domnames[idx] for idx in pair]
        pair_dom_dbids = dom_dbids[pair,] # [dom_dbids[idx] for idx in pair]

        chn1 = chain_dom_structs[pair_domnames[0]]
        chn2 = chain_dom_structs[pair_domnames[1]]

        if domains_in_contact(chn1, chn2):

            order = np.argsort(pair_dom_dbids)

            # the domain in the pair that will be used for superposition
            ref_dom = pair_domnames[order[0]]
            ref_dom_str = deepcopy(chain_dom_structs[ref_dom])

            # the tag-along
            tag_along_dom = pair_domnames[order[1]]
            tag_along_dom_str = deepcopy(chain_dom_structs[tag_along_dom])

            dompair_id = sep.join(pair_dom_dbids[order])
            domnames_in_order = domname_sep.join(pair_domnames[order])

            if dompair_id not in result.keys():
                result[dompair_id] = dict()

                # this is the 'master' domain to which all instances of the same dbid pair must be aligned,
                #  aligning on the domain with domid of lower lexical value

                # result[dompair_id]['ref_dom_str'] = ref_dom_str # no need to store the structure itself, only the id
                result[dompair_id]['ref_domid'] = ref_dom

                com_reference = ref_dom_str.center_of_mass()
                com_tag_along = tag_along_dom_str.center_of_mass()

                v = com_tag_along - com_reference

                result[dompair_id]['vectors'] = v.reshape((1,3))

                # these _can_ be omitted if space becomes an issue
                # result[dompair_id]['ref_domain_pair'] = domname_sep.join((ref_dom, tag_along_dom))

                result[dompair_id]['aligned_domain_pairs'] = list()
                result[dompair_id]['aligned_domain_pairs'].append(domnames_in_order)

            else:
                # align on ref_domain
                ref_dom_fname = os.path.join(dom_pdb_dir, ref_dom)
                master_ref_dom_fname = os.path.join(dom_pdb_dir, result[dompair_id]['ref_domid'])

                # TM-align to get rotation+translation matrix+vector
                tmscore, aln, transl, rot, errflag = tmalign.tmalign(ref_dom_fname, master_ref_dom_fname, args=tmalign_args, exe=tmalign_bin)

                if errflag > 0:
                    print("Error in tmalign step; skipping pair", pair_domnames, dompair_id)
                    continue

                if warn and tmscore <= tmscore_warn_thresh:
                    print("WARNING: TM-align("+ref_dom_str.id+", "+result[dompair_id]['ref_domid']+") = "+str(tmscore)+" !", file=sys.stderr)
                #transform coords of current domain pair
                for a in ref_dom_str.get_atoms():
                    a.transform(np.transpose(rot), transl)
                for a in tag_along_dom_str.get_atoms():
                    a.transform(np.transpose(rot), transl)

                # compute vectors and other data
                com_reference = ref_dom_str.center_of_mass()
                com_tag_along = tag_along_dom_str.center_of_mass()

                v = com_tag_along - com_reference
                result[dompair_id]['vectors'] = np.concatenate( [ result[dompair_id]['vectors'], v.reshape(1,3) ] )

                result[dompair_id]['aligned_domain_pairs'].append(domnames_in_order)

    for f,d in zip(dom_pdb_fnames, domnames):
        f_in_tmp = os.path.join(tmpdir, d)
        os.remove(f_in_tmp)

pickle.dump(result, open(iv_out_f, 'wb'))
print('All completed')

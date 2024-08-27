#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 11:02:38 2023

Main domain pair evaluation script for TED domains

@author: Shaun M Kandathil
"""
###############################################################################
## IMPORTS

import sys,os
import shutil
import pickle
import subprocess
from itertools import combinations, chain
from copy import deepcopy

import numpy as np
from Bio.PDB import PDBParser, NeighborSearch

import tmalign

###############################################################################
## GLOBALS

basedir = os.path.join(os.environ['HOME'], 'afdb_domain', 'tools', 'dompair-eval-afdb')

# I say tgz but actually is zip
tgz_dir = os.path.join(os.environ['HOME'], 'afdb_domain', 'afdb_domain_ext', 'results', 'proteome_consensus_domains_zip')
tgz_prefix = "proteome-tax_id-"
tgz_suffix = "_v4.consensus_domains" # .tar.gz/.zip added elsewhere
# proteome-tax_id-53326-6_v4.consensus_domains.zip

refdom_d_fname = "/cluster/project9/afdb_domain_ext/results/dompair-eval-inputs/ref_domains_by_tedid.pkl" # sys.argv[2] #
refdom_pdb_dir = "/cluster/project9/afdb_domain_ext/results/dompair-eval-inputs/master-ref-domain-pdbs"

tmpdir = '/scratch0/' + os.environ['LOGNAME']
os.makedirs(tmpdir, exist_ok=True)

out_dir = "/cluster/project9/afdb_domain_ext/results/dompair-eval-outputs-v2"
os.makedirs(out_dir, exist_ok=True)

# atom distance thresholds in Angstrom for Cbeta and any-heavy-atom contacts
contact_params = {'cbeta': 8.0, 'heavy': 3.0}

# contact definition based on C-beta atom ('cbeta') or any heavy atom ('heavy')
cont_type = 'cbeta'

# min number of contacts between 2 domains to say that they are in contact
min_contacts = 3

# separator used to join dbids and domnames
sep = '-'
domname_sep = ':'

# generate domain pair dict? Adds runtime and will take space+inodes when dumped
gen_dompair_info = False

# path to tmalign executable
tmalign_bin = os.path.join(os.environ['HOME'], "software_staging", "tmalign", "tmalign")
if os.path.exists(tmpdir):
    shutil.copy2(tmalign_bin, tmpdir+"/tmalign")
    tmalign_bin = tmpdir+"/tmalign"

#chaindomdic (lookup-dict) prefix and suffix
cdd_prefix = os.path.join("/cluster/project9/afdb_domain_ext/results/dompair-eval-inputs-chunks",'lookup-dict-chunk')
cdd_suffix = '.pkl'

# string in pdb filenames
pdb_mid_str = '-model_v4_'

# output filename base
iv_out_f_base = '/'.join([out_dir,'dompair_interaction_vectors_chunk'])

# tmalign args (-m is always specified. NEVER set -outfmt)
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

# os.chdir('/cluster/project9/afdb_domain_ext/results/dompair-eval-inputs')
chunk_id = sys.argv[1] #"000001" # or just '1' #
chunk_id = chunk_id.zfill(6)
cdd_fname = cdd_prefix + chunk_id + cdd_suffix

iv_out_f = iv_out_f_base + chunk_id + ".pkl"

print('processing', cdd_fname)

if not os.path.exists(refdom_pdb_dir):
    print('ERROR: refdom_pdb_dir', refdom_pdb_dir, "not accessible/doesn't exist..", file=sys.stderr)
    sys.exit(1)

chaindomdic = pickle.load(open(cdd_fname,'rb'))
refdom_d = pickle.load(open(refdom_d_fname, 'rb'))
# '3.10.450.140': {'taxid': '562-28', 'chain': 'AF-A0A376W4B5-F1-model_v4_TED02', 'plddt': '79.3871'}

print("chaindomdic and refdom_d successfully loaded")
dbids = 'tedids'

# first filter the chaindomdic to only retain chains with 2 or more domains.
# SMK Andy has already filtered these out.
#chaindomdic = {k:v for k,v in chaindomdic.items() if len(v["domids"]) > 1}

proteomelist = chaindomdic.keys()

result = dict()
chaincount=0

parser = PDBParser(QUIET=True)

for k1 in proteomelist:
    # copy the entire tgz to tmpdir (scratch)
    tgz_base     = tgz_prefix + k1 + tgz_suffix
    tgz_path     = tgz_dir + '/' + tgz_base + '.zip'
    tgz_f_in_tmp = tmpdir  + "/" + tgz_base + '.zip'
    if not os.path.exists(tgz_f_in_tmp):
        shutil.copy(tgz_path, tgz_f_in_tmp)

    for k in chaindomdic[k1].keys(): # k is a chain id
        chaincount+=1
        chain_dom_structs = dict() # name-addressable list

        try:
            domnames  = chaindomdic[k1][k]['domids']
            dom_dbids = chaindomdic[k1][k][dbids]
            dom_pdb_dir = tmpdir # os.path.join(tmpdir, tgz_base)
            dom_pdb_fnames = [k+pdb_mid_str+a+'.pdb' for a in domnames]
            domnames = np.asarray(dom_pdb_fnames)
        except Exception as e:
            print(e)
            print(chaindomdic[k1][k])
            sys.exit(1)


        ## tar will output files like proteome-tax_id-8010-4_v4.consensus_domains/AF-A0A6Q2XI59-F1-model_v4_TED01.pdb

        # extract all domains for the current chain
        # again, it says tar but we actually use unzip
        #os.chdir(tmpdir)
        # tarcmd = ['tar', '-C', tmpdir, '-xf', tgz_f_in_tmp] + dom_pdb_fnames
        tarcmd = ['unzip', '-q', '-o', '-d', tmpdir, tgz_f_in_tmp] + dom_pdb_fnames
        tarproc = subprocess.run(args=tarcmd, capture_output=False)
        if tarproc.returncode != 0:
            print('unzip error:', file=sys.stderr)
            # print('in dir '+tmpdir, file=sys.stderr)
            print('command:',tarcmd, file=sys.stderr)

            sys.exit(1)
        os.chdir(tmpdir)
        for f,d in zip(dom_pdb_fnames, domnames):
            # read
            chain_dom_structs[d] = parser.get_structure(d, f)

        os.chdir(basedir)

        ndom = len(domnames)

        # find contacting pairs of domains
        combns = combinations(list(range(ndom)), 2)

        for pair in combns:
            # i,j = pair # indices into each of the 2 arrays in chaindomdic[k1][k]

            pair_domnames  =  domnames[pair,]  # [domnames[idx] for idx in pair]
            pair_dom_dbids = dom_dbids[pair,] # [dom_dbids[idx] for idx in pair]

            chn1 = chain_dom_structs[pair_domnames[0]]
            chn2 = chain_dom_structs[pair_domnames[1]]

            if domains_in_contact(chn1, chn2):

                order = np.argsort(pair_dom_dbids)

                # the domain in the current pair that will be used for superposition
                ref_dom = pair_domnames[order[0]]
                ref_domid = pair_dom_dbids[order[0]]
                ref_dom_str = deepcopy(chain_dom_structs[ref_dom])

                # the tag-along in the current pair
                tag_along_dom = pair_domnames[order[1]]
                tag_along_dom_str = deepcopy(chain_dom_structs[tag_along_dom])

                # the pre-determined 'master reference' domain for the dbid
                # fetch from refdom_d
                master_ref_dom_fname = refdom_d[ref_domid]['chain']
                master_ref_dom_taxid = refdom_d[ref_domid]['taxid']
                # tgz_base     = tgz_prefix + master_ref_dom_taxid + tgz_suffix
                master_ref_dom_pdb_path = '/'.join([refdom_pdb_dir, master_ref_dom_fname+'.pdb'])

                dompair_id = sep.join(pair_dom_dbids[order])
                domnames_in_order = domname_sep.join(pair_domnames[order]).replace('.pdb','')

                if not dompair_id in result.keys():
                    result[dompair_id] = dict()

                    # # these _can_ be omitted if space becomes an issue
                    result[dompair_id]['ref_domain_pair'] = domname_sep.join((ref_dom, tag_along_dom)).replace('.pdb','')

                    result[dompair_id]['aligned_domain_pairs'] = list()

                # align on ref_domain
                ref_dom_fname = os.path.join(tmpdir, ref_dom)

                # TM-align to get rotation+translation matrix+vector
                tmscore, aln, transl, rot, errflag = tmalign.tmalign(ref_dom_fname, master_ref_dom_pdb_path, args=tmalign_args, exe=tmalign_bin)

                if errflag > 0:
                    print("Error in tmalign step; skipping pair", ref_dom_fname, master_ref_dom_pdb_path, dompair_id)
                    continue

                if warn and tmscore <= tmscore_warn_thresh:
                    print("WARNING: TM-align("+ref_dom_str.id+", "+ master_ref_dom_fname +") = "+str(tmscore)+" !", file=sys.stderr)

                # transform coords of current domain pair
                for a in ref_dom_str.get_atoms():
                    a.transform(np.transpose(rot), transl)
                for a in tag_along_dom_str.get_atoms():
                    a.transform(np.transpose(rot), transl)

                # compute vectors and other data
                com_reference = ref_dom_str.center_of_mass()
                com_tag_along = tag_along_dom_str.center_of_mass()

                v = com_tag_along - com_reference
                # save the raw vectors - conversion to unit vecs etc happens downstream.
                if 'vectors' in result[dompair_id].keys():
                    result[dompair_id]['vectors'] = np.concatenate( [ result[dompair_id]['vectors'], v.reshape(1,3) ] )
                else:
                    result[dompair_id]['vectors'] = v.reshape((1,3))

                result[dompair_id]['aligned_domain_pairs'].append(domnames_in_order)

        for f,d in zip(dom_pdb_fnames, domnames):
            f_in_tmp = os.path.join(tmpdir, f)
            os.remove(f_in_tmp)

        # if chaincount % 10 == 0:
        #     print(chaincount)

    os.remove(tgz_f_in_tmp)


pickle.dump(result, open(iv_out_f, 'wb'))
print('All completed')

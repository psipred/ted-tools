#!/usr/bin/env python3

# Utilities for interdomain PAE data 
# Shaun Kandathil, 2023

import urllib.request
import json
import numpy as np

pae_url_base = "https://alphafold.ebi.ac.uk/files/"

def chopping_to_res(chopping, range_sep='-', discont_sep='_'):
    '''
    Take a TED chopping string and return the complete list of residue indices implied by it.

    '''
    subranges = chopping.split(discont_sep)
    result = list()
    for subrange in subranges:
        start, end = subrange.split(range_sep)
        result.extend(list(range(int(start), int(end)+1)))

    return np.asarray(result)

def gen_np_index(align_res, score_res, both=True, minus1=True):
    '''
    Return NumPy indices tuple for a pair of domains, for use with pAE matrices.
    
    @param align_res: list of ints; residue numbers for the 'aligned' residues, usually taken from one domain
    @param score_res: list of ints; residue numbers for the 'scored' residues, usually taken from one domain
    @param both: boolean; if True, the output tuple has indices for both triangles of the pAE matrix.
    @param minus1: boolean; if True, `align_res` and `score_res` are treated as 1-indexed residue numbers, and 1 is subtracted from them to yield the 0-indexed version in the output.
    '''
    
    lar = len(align_res)
    lsr = len(score_res)

    if minus1:
        align_res -= 1
        score_res -= 1
    # lower triangle
    align_idx1 = np.repeat(align_res, lsr)
    score_idx1 = np.tile(score_res, lar)

    if not both:
        return (align_idx1, score_idx1)
    # upper triangle
    align_idx2 = np.repeat(score_res, lar)
    score_idx2 = np.tile(align_res, lsr)
    
    return (np.append(align_idx1, align_idx2), np.append(score_idx1, score_idx2))


def interdom_pae_summary(pae, align_res, score_res, both=True, minus1=True, func=np.median):
    '''
    Extract pAE submatrix, and return the result of applying `func` to it.

    @param pae: np.array of ints with shape [N, N]; the pAE matrix.
    @param align_res: see `gen_np_index`
    @param score_res: see `gen_np_index`
    @param both: see `gen_np_index`
    @param minus1: see `gen_np_index`
    @param func: callable; a function to be applied to the extracted pAE submatrix. Defaults to np.mean.
    '''
    idx = gen_np_index(align_res, score_res, both, minus1)
    return(func(pae[idx]))

def get_pae_from_afdb_id(id_, url_base=pae_url_base):
    paeurl = url_base + id_ + '-predicted_aligned_error_v4.json'
    with urllib.request.urlopen(paeurl) as url:
        pae = np.asarray(json.load(url)[0]['predicted_aligned_error'])

    return pae

if __name__ == '__main__':
    
    afdb_id = "AF-A0A086H0Z9-F1"
    chopping_d1 = "3-153_165-226"
    chopping_d2 = "262-397"

    res_d1 = chopping_to_res(chopping_d1)
    res_d2 = chopping_to_res(chopping_d2)
    
    pae = get_pae_from_afdb_id(afdb_id)
    
    score_symm = interdom_pae_summary(pae, res_d1, res_d2)
    score_d1d2 = interdom_pae_summary(pae, res_d1, res_d2, both=False)
    score_d2d1 = interdom_pae_summary(pae, res_d2, res_d1, both=False)
    print(afdb_id, chopping_d1, chopping_d2)
    print("D1:", chopping_d1)
    print("D2:", chopping_d2)
    print("Symmetric score:", score_symm)
    print("A(D1) S(D2) score:", score_d1d2)
    print("A(D2) S(D1) score:", score_d2d1)

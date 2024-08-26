from itertools import groupby
from operator import itemgetter
from natsort import natsorted

import numpy as np
import pandas as pd

# This file is a part of TED: The Encyclopedia of Domains. If you utilize or reference any content from this file, 
# please cite the following paper:

# Lau et al., 2024. Exploring structural diversity across the protein universe with The Encyclopedia of Domains.

def read_chopping(file: str, headings: list):
    """ Reads a chopping file into a pandas dataframe.
    
    Intended chopping format:
    One target per row, tab separated values, fields: [target, md5, nres, ndom, chopping(doms delimited by , discontinuous delimited by _ ranges by -), scores (either single float or csv)]
    
    AF-O94910-F1-model_v4   01b0f479da155fad9df6f696723449e9        1474    6       27-137,144-396,474-533,534-639,640-838,847-1131 1028.89,81.89,1017.89,1035.82,1053.72,1075.68
    AF-O94779-F1-model_v4   0ed866265ad5f1bf2496dea7dba8af24        1100    10      95-192,193-296,297-387,388-477,478-569,570-671,672-768,769-870,871-972,973-1064 1025.21,1023.13,1022.06,1022.41,1018.82,1009.06,1020.49,1014,1025.89,1011.89

    Ignores the 'Done.' entry at the bottom of chopping files if its present.

    Args:
        file        (str)               Path to chopping file.
        headings    (list)              List of heading names to use.
    
    Returns:
        out         (dataframe)             Pandas dataframe object.
    
    """
    df = pd.read_csv(file, delimiter='\t', header=None, names=headings)
    df = df[~df['target'].str.match('Done.')]
    return df

def domstr_to_ranges(chopping: str, delim_dom: str=',', delim_discon: str='_', delim_range: str='-', offset: int=0):
    """ Parses a domain string into a list of residue ranges.

    Args:
        chopping        (str)               Domain string, e.g. '1-100_200-250,251-300'
        delim_dom       (str, optional)     Delimiter for domains. Defaults to ','.
        delim_discon    (str, optional)     Delimiter for discontinuous segments. Defaults to '_'.
        delim_range     (str, optional)     Delimiter for residue ranges. Defaults to '-'.

    Returns:
        out     (list)  List of domain ranges in format:
                        [[domain_id, start, end]], e.g.
                        [[1, 71, 83], [1, 143, 227], [2, 84, 142], [2, 228, 244]]
    """
    
    domains = []
    
    dom_count = 1
    for _, dom in enumerate(chopping.split(delim_dom)):
        doms = [d for d in dom.split(delim_discon)]
        
        for d in doms:
            if '-' in d:
                start, end = d.split(delim_range)
                start, end = int(start), int(end)
                
                start += offset
                end += offset
                
                # Below start shifts by -1 so that residue 1 is index 0. 
                # assert start > 0, f"Start residue cannot be less than 1. Got residue {start}. Note: You can use the --offset_f1 and --offset_f2 flags to apply an offset to the choppings for each file."
                # assert end > start, f"End residue cannot be less than 1 and lower than start residue. Got residue {end}. Note: You can use the --offset_f1 and --offset_f2 flags to apply an offset to the choppings for each file."

                domains.append([dom_count, start, end])
                
        dom_count += 1

    return domains

def domstr_to_assignment_by_resi(chopping: str, resi: np.ndarray, offset: int = 0):
    """ Convert a domstr (domain string) into an array of domain indices

    Args:
        chopping        (str)               Domain string, e.g. '1-100_200-250,251-300'
                                            Residues expected to start from 1. 
        resi            (np.ndarray)        NumPy array of target residues. 
        offset          (int)               Optional, shift the start and end residues by n.

    Returns:
                        (np.ndarray)         NumPy array representing the domain index of each residue.
                                            Ranges from 0 (NDR) to ndom.    
    """

    assignment = np.zeros(len(resi))

    if chopping != '0' and chopping != 'NULL':
        domains = domstr_to_ranges(chopping, offset=offset)

        for d in domains:
            dom_id, start, end = d

            # Below start shifts by -1 so that residue 1 is index 0. 
            # assert start >= 0, f"Start index cannot be less than 0. Got residue {start}."
            # assert end > start, f"End residue cannot be less than the start index. Got residue {end}."

            idx = np.where((resi >= start) & (resi <= end))[0]
            assignment[idx] = dom_id

    return assignment.astype(np.int64)

def domstr_to_assignment(chopping: str, nres: int, offset: int = 0):
    """ Convert a domstr (domain string) into an array of domain indices

    Args:
        chopping        (str)               Domain string, e.g. '1-100_200-250,251-300'
                                            Residues expected to start from 1. 
        nres            (int)               The length of the target, N. 
        offset          (int)               Optional, shift the start and end residues by n.

    Returns:
                        (np.ndarray)         NumPy array representing the domain index of each residue.
                                            Ranges from 0 (NDR) to ndom.    
    """
    
    assignment = np.zeros(nres)

    if chopping != '0' and chopping != 'NULL':
        domains = domstr_to_ranges(chopping, offset=offset)

        for d in domains:
            dom_id, start, end = d
            
            # Below start shifts by -1 so that residue 1 is index 0. 
            # assert start >= 0, f"Start index cannot be less than 0. Got residue {start}."
            # assert end > start, f"End residue cannot be less than the start index. Got residue {end}."

            assignment[start-1:end] = dom_id

    return assignment.astype(np.int64)

def assignment_to_domstr(assign: np.ndarray, ri: np.ndarray, delim_dom: str=',', delim_discon: str='_', delim_range: str='-'):

    unique_ids = np.unique(assign[assign != 0])
    ndom = len(unique_ids)

    dom_str = []
    for d in unique_ids:
        dom_resi = ri[assign == d]

        strs = []
        for k, g in groupby(enumerate(dom_resi), lambda ix: ix[0] - ix[1]):
            consecs = list(map(itemgetter(1), g))
            if len(consecs) > 1:
                strs.append("{}{}{}".format(str(consecs[0]), delim_range, str(consecs[-1])))
            else:
                strs.append(str(consecs[0]))

        dom_str.append(delim_discon.join(strs))

    new_domstr = delim_dom.join(natsorted(dom_str))
    if new_domstr == '':
        new_domstr = '0'

    return new_domstr, ndom

def tabulate_tensor(x, reverse=False):
    if reverse:
        return {x.bincount()[n].item(): n.item() for n in x.unique()}
    else:
        return {n.item(): x.bincount()[n].item() for n in x.unique()}

def filter_domstr(domstr, segment_threshold=5, domain_threshold=25):
    new_segments = []
    domain_nres = 0
    for segment in domstr.split('_'):
        if '-' in segment:
            start, end = segment.split('-')
            if (int(end) - int(start) + 1) > segment_threshold:
                new_segments.append(segment)
                domain_nres += int(end) - int(start) + 1
    
    if domain_nres >= domain_threshold and len(new_segments) > 0:                
        domstr = '_'.join(natsorted(new_segments))
    else:
        domstr = '0'
        
    return domstr
import numpy as np
import networkx as nx

from itertools import product
from string import ascii_lowercase
from natsort import natsorted

from score_utils import domstr_to_assignment, assignment_to_domstr, domstr_to_ranges, filter_domstr

# This file is a part of TED: The Encyclopedia of Domains. If you utilize or reference any content from this file, 
# please cite the following paper:

# Lau et al., 2024. Exploring structural diversity across the protein universe with The Encyclopedia of Domains.

# Replaced aa so that it handles longer than 26 domains
aa = [''.join(combo) for x in range(1, 6) for combo in product(ascii_lowercase, repeat=x)]
# range(1,4) handles up to 18,278 
# range(1,6) handles up to 12,356,630 domains

def calculate_domain_consensus(choppings: list, nres: int, iou_threshold: float=0.7, consensus_levels: list=[4,3,2]):

    ri = np.arange(1, nres+1)

    assignments = []
    assign_by_str = {}
    unique_doms = []
    for i, c in enumerate(choppings):
        assign = domstr_to_assignment(c, nres, offset=1)
        labels = np.char.mod(f'{aa[i]}%d', assign)
        doms = np.unique(labels)

        assignments.append(assign)
        assign_by_str[aa[i]] = labels
        unique_doms.append(doms)
        
    high, medium, low = [], [], []
    high_domstr, med_domstr, low_domstr = [], [], []

    if len(unique_doms) > 0:
        unique_doms = np.hstack(unique_doms)
        
        unique_doms = np.delete(
            unique_doms, np.where(np.char.find(np.char.lower(unique_doms), '0') > -1)[0]
        )

        # Above generates three outputs:
        # assignments - this is just [0,0,0,0,1,1,1,1,2,2,2, ...] for each method
        # assign_by_str - this is a dictionary containing:
        #   {'a': [a0,a0,a0,a0,a1,a1,a1,a1,a2,a2,a2, ...]
        #    'b': [b0,b0,b0,b0,b1,b1,b1,b2,b2,b2,b2, ...] ... }
        #   Where each residue's assignment is prepended with a letter unique to that method
        #   The assignment of a method can be accessed using assignment_by_str[a], etc. 
        # unique_doms - this is a list of unique non-NDR domains, e.g. a1-9, no a0. 

        ndoms = unique_doms.shape[0]

        # Calculate consensus matrix:

        # Generate an N-by-N adjacency matrix and iterate over pairs of domains and calculate IoU
        assign_by_id = {}

        adj_dom = np.zeros((ndoms, ndoms))
        for i in range(ndoms):
            di = unique_doms[i]
            ci = di[0]
            di_assign = assign_by_str[ci] == di
            assign_by_id[di] = di_assign
            
            for j in range(ndoms):
                dj = unique_doms[j]
                cj = dj[0]
                
                dj_assign = assign_by_str[cj] == dj
                
                intersect = np.sum(np.logical_and(di_assign, dj_assign))
                union = np.sum(np.logical_or(di_assign, dj_assign))
                iou = intersect / union
                
                adj_dom[i, j] = iou
            
        # For each domain, sum up the number of corresponding domains meeting an IoU threshold
        # n_consensus = np.sum(adj_dom >= iou_threshold, axis=1)
        # for d, n in zip(unique_doms, n_consensus):
        #     print(d, n)
      
        G = (adj_dom >= iou_threshold).astype(np.int32)
        components = list(nx.connected_components(nx.from_numpy_array(G)))
        for i in components:
            
            if len(i) > 0:
                doms = [assign_by_id[unique_doms[idx]] for idx in i]
                consensus = np.sum(np.vstack(doms), axis=0)

            # Set non-max values to 0 so that we get the intersect of the best overlapping parts
            if np.max(consensus) > 0:
                consensus[consensus != np.max(consensus)] = 0
                
            if len(i) >= consensus_levels[0]:
                domstr = assignment_to_domstr(
                    (consensus >= consensus_levels[0]).astype(np.int32), ri)[0]

                if domstr != '0':
                    high.extend(domstr_to_ranges(domstr))
                    high_domstr.append(domstr)
            
            if len(i) >= consensus_levels[1] and len(i) < consensus_levels[0]:
                domstr = assignment_to_domstr(
                    np.where((consensus >= consensus_levels[1]) & (consensus < consensus_levels[0]), 1, 0), ri)[0]

                if domstr != '0':
                    medium.extend(domstr_to_ranges(domstr))
                    med_domstr.append(domstr)
                
            if len(i) >= consensus_levels[2] and len(i) < consensus_levels[1]:
                domstr = assignment_to_domstr(
                    np.where((consensus >= consensus_levels[2]) & (consensus < consensus_levels[1]), 1, 0), ri)[0]
                
                if domstr != '0':
                    low.extend(domstr_to_ranges(domstr))
                    low_domstr.append(domstr)
    
    if len(high_domstr) > 0:
        n_high = len(high_domstr)
        high_domstr = ','.join(natsorted(high_domstr))
    else:
        high_domstr = 'na'
        n_high = 0
        
    if len(med_domstr) > 0:
        n_med = len(med_domstr)
        med_domstr = ','.join(natsorted(med_domstr))
    else:
        med_domstr = 'na'
        n_med = 0
        
    if len(low_domstr) > 0:
        n_low = len(low_domstr)
        low_domstr = ','.join(natsorted(low_domstr))
    else:
        low_domstr = 'na'
        n_low = 0

    return (high, medium, low), (high_domstr, med_domstr, low_domstr), (n_high, n_med, n_low)


def calculate_domain_consensus_redundant(choppings: list, nres: int, iou_threshold: float=0.7, consensus_levels: list=[4,3,2]):
    # for sequence-redundant consensus

    ri = np.arange(1, nres+1)

    assignments = []
    assign_by_str = {}
    unique_doms = []
    for i, c in enumerate(choppings):
        assign = domstr_to_assignment(c, nres)
        labels = np.char.mod(f'{aa[i]}-%d', assign)
        doms = np.unique(labels)

        assignments.append(assign)
        assign_by_str[aa[i]] = labels
        unique_doms.append(doms)

    high, medium, low = [], [], []
    high_domstr, med_domstr, low_domstr = [], [], []
    high_assign = []

    if len(unique_doms) > 0:
        unique_doms = np.hstack(unique_doms)
        
        unique_doms = np.delete(
            unique_doms, np.where(np.char.find(np.char.lower(unique_doms), '0') > -1)[0]
        )

        # Above generates three outputs:
        # assignments - this is just [0,0,0,0,1,1,1,1,2,2,2, ...] for each method
        # assign_by_str - this is a dictionary containing:
        #   {'a': [a0,a0,a0,a0,a1,a1,a1,a1,a2,a2,a2, ...]
        #    'b': [b0,b0,b0,b0,b1,b1,b1,b2,b2,b2,b2, ...] ... }
        #   Where each residue's assignment is prepended with a letter unique to that method
        #   The assignment of a method can be accessed using assignment_by_str[a], etc. 
        # unique_doms - this is a list of unique non-NDR domains, e.g. a1-9, no a0. 

        ndoms = unique_doms.shape[0]

        # Calculate consensus matrix:

        # Generate an N-by-N adjacency matrix and iterate over pairs of domains and calculate IoU
        assign_by_id = {}

        adj_dom = np.zeros((ndoms, ndoms))
        for i in range(ndoms):
            di = unique_doms[i]
            ci = di.split('-')[0]
            di_assign = assign_by_str[ci] == di
            assign_by_id[di] = di_assign

            for j in range(ndoms):
                dj = unique_doms[j]
                cj = dj.split('-')[0]
                
                dj_assign = assign_by_str[cj] == dj
                
                intersect = np.sum(np.logical_and(di_assign, dj_assign))
                union = np.sum(np.logical_or(di_assign, dj_assign))
            
                if union == 0 or np.isnan(union):
                    iou = 0
                else:
                    iou = intersect / union
                
                adj_dom[i, j] = iou
      
        G = (adj_dom >= iou_threshold).astype(np.int32)
        # components = list(nx.connected_components(nx.from_numpy_array(G)))
        components = sorted(nx.connected_components(nx.from_numpy_array(G)), key=len, reverse=True)
        
        high_comps = [c for c in components if len(c) >= consensus_levels[0]]
        med_comps = [c for c in components if len(c) >= consensus_levels[1] and len(c) < consensus_levels[0]]
        low_comps = [c for c in components if len(c) >= consensus_levels[2] and len(c) < consensus_levels[1]]

        # Iterate over all high consensus domains first then medium. This allows us to check that 
        # medium consensus domains do not overlap with the high consensus ranges.
        
        if len(high_comps) > 0:
            for i in high_comps:
                doms = [assign_by_id[unique_doms[idx]] for idx in i]
                consensus = np.sum(np.vstack(doms), axis=0)
                
                # Set non-max values to 0 so that we get the intersect of the best overlapping parts
                consensus[consensus != np.max(consensus)] = 0
                
                # Convert into binary index to get intersect of all assignments in this component
                assignment = (consensus >= consensus_levels[0]).astype(np.int32)
                domstr = filter_domstr(assignment_to_domstr(assignment, ri)[0])
                
                # Hack: only save the assignment if its not '0'
                if domstr != '0':
                    high.extend(domstr_to_ranges(domstr))
                    high_domstr.append(domstr)
                    high_assign.append(assignment)
                    
            # Get the overall high consensus assignment as a binary array
            # 1 where a residue is in a high consensus domain, 0 otherwise
            if len(high_assign) > 1:
                high_assign = (np.sum(np.vstack(high_assign), axis=0) > 0).astype(int)
            else:
                high_assign = assignment
        
        if len(med_comps) > 0:
            for i in med_comps:
                doms = [assign_by_id[unique_doms[idx]] for idx in i]
                consensus = np.sum(np.vstack(doms), axis=0)
                
                # Set non-max values to 0 so that we get the intersect of the best overlapping parts
                consensus[consensus != np.max(consensus)] = 0
                
                # Convert into binary index to get intersect of all assignments in this component
                assignment = np.where((consensus >= consensus_levels[1]) & (consensus < consensus_levels[0]), 1, 0)

                # Remove any domains that overlap with high consensus ones
                if len(high_assign) > 0:
                    if np.sum(assignment * (high_assign == 0)) > 0:
                        continue
                    
                # convert to domain string format
                domstr = filter_domstr(assignment_to_domstr(assignment, ri)[0])

                # Hack: only save the assignment if its not '0'
                if domstr != '0':
                    medium.extend(domstr_to_ranges(domstr))
                    med_domstr.append(domstr)

        if len(low_comps) > 0:
            for i in low_comps:
                doms = [assign_by_id[unique_doms[idx]] for idx in i]
                consensus = np.sum(np.vstack(doms), axis=0)
                
                # Set non-max values to 0 so that we get the intersect of the best overlapping parts
                consensus[consensus != np.max(consensus)] = 0
                
                # Convert into binary index to get intersect of all assignments in this component
                assignment = np.where((consensus >= consensus_levels[2]) & (consensus < consensus_levels[1]), 1, 0)

                # no filtering for low consensus domains as we don't care about overlap here
                # # Remove any residues in high consensus domains
                # if len(high_assign) > 0:
                #     assignment = assignment * (high_assign == 0)
                
                # convert to domain string format
                domstr = filter_domstr(assignment_to_domstr(assignment, ri)[0])

                # Hack: only save the assignment if its not '0'
                if domstr != '0':
                    low.extend(domstr_to_ranges(domstr))
                    low_domstr.append(domstr)
    
    if len(high_domstr) > 0:
        n_high = len(high_domstr)
        high_domstr = ','.join(natsorted(high_domstr))
    else:
        high_domstr = 'na'
        n_high = 0
        
    if len(med_domstr) > 0:
        n_med = len(med_domstr)
        med_domstr = ','.join(natsorted(med_domstr))
    else:
        med_domstr = 'na'
        n_med = 0
        
    if len(low_domstr) > 0:
        n_low = len(low_domstr)
        low_domstr = ','.join(natsorted(low_domstr))
    else:
        low_domstr = 'na'
        n_low = 0

    return (high, medium, low), (high_domstr, med_domstr, low_domstr), (n_high, n_med, n_low)


if __name__ == "__main__":
    
    # For testing
    choppings = [
        '53-145,210-310_377-384,311-376_385-398,399-489,490-571,572-746,747-918,919-1121',
        '54-138,166-196,210-304,306-394,409-483,490-569,576-746,750-914,920-1118',
        '53-145,223-305,306-399,400-485,486-573,574-749,750-918,919-1121',
        ]

    nres = 1184
    
    _, domstrs, counts = calculate_domain_consensus(choppings, nres, iou_threshold=0.7, consensus_levels=[3,2,1])
    
    high, medium, low = domstrs
    nhigh, nmedium, nlow = counts
    
    print("Nres: ", nres)
    print("Low consensus: ", low)
    print("Medium consensus: ", medium)
    print("High consensus: ", high)
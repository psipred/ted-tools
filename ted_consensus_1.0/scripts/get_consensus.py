import os
import time
import argparse

import numpy as np

from score_utils import domstr_to_ranges
from domain_consensus import calculate_domain_consensus

# This file is a part of TED: The Encyclopedia of Domains. If you utilize or reference any content from this file, 
# please cite the following paper:

# Lau et al., 2024. Exploring structural diversity across the protein universe with The Encyclopedia of Domains.

scriptdir = os.path.dirname(os.path.realpath(__file__))

def read_chopping(file):
    with open(file, 'r') as f:
        contents = []
        targets = []
        for line in f:
            if line[:5] != 'Done.':
                line = line.rstrip('\n').split('\t')
                
                assert len(line) == 6, f"Expected 6 fields in chopping file {file} in format: 'target,md5,nres,ndom,chopping,score'"
                
                targets.append(line[0])
                contents.append(line)
        
    return targets, np.array(contents)

def main():
    # Read the config file
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--choppings", type=str, nargs="+", required=True, help="Pass a list of files containing domain choppings.")
    parser.add_argument("-o", "--output", type=str, required=True, help="File to save consensus domains to.")
    parser.add_argument("--overlap", type=float, required=False, default=0.7, help="IoU threshold to use for determining matches.")
    args = parser.parse_args()
    
    start_time = time.time()
    
    # Open the chopping files
    targets = []
    chopping_files = {}
    for f in args.choppings: 
        start_read = time.time()
        bn, _ = os.path.splitext(os.path.basename(f))
        chopping_targets, chopping_files[bn] = read_chopping(f)
        targets.extend(chopping_targets)

        print("Read file {} in {:.2f} seconds".format(bn, time.time() - start_read))
        
    print("Finished file importing in {:.2f} seconds".format(time.time() - start_time))
        
    # Get unique list of targets from all chopping files
    unique_targets = list(set(targets))
    
    with open(args.output, 'w+') as f:
        for _, target in enumerate(unique_targets):
            
            dom_method = []
            method_choppings = []
            for k, v in chopping_files.items():
                match = v[v[:,0] == target]
                
                # Check that only one match can be made
                assert len(match) > 0, f"Entry not found for {target} in chopping file {k}"
                assert len(match) < 2, f"Multiple matches for target {target} in chopping file {k}."
                
                # Format domain string into a pymol command
                target, md5, nres, ndom, chopping, score = match[0]
                nres, ndom, score = int(nres), int(ndom), float(score)

                if chopping not in ['0','NULL','NO_SS']:
                    dom_method.append(domstr_to_ranges(chopping))
                    method_choppings.append(chopping)
                else:
                    dom_method.append([[0, 1, nres-1]])

            # Calculate the consensus between the multiple assignments:
            _, consensus_domstr, consensus_counts = calculate_domain_consensus(method_choppings, nres, consensus_levels=[3,2,1], iou_threshold=args.overlap)
            high, medium, low = consensus_domstr
            nhigh, nmed, nlow = consensus_counts
            
            f.write("{}\t{}\t{:.0f}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                target, md5, nres, nhigh, nmed, nlow, high, medium, low,
            ))


if __name__ == "__main__":
    main()
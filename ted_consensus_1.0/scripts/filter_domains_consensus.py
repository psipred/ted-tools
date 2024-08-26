import argparse
import os
from natsort import natsorted

# This file is a part of TED: The Encyclopedia of Domains. If you utilize or reference any content from this file, 
# please cite the following paper:

# Lau et al., 2024. Exploring structural diversity across the protein universe with The Encyclopedia of Domains.

def filter_domains(domains, args):
                
    new_ndom = 0
    new_dom_chopping = []
    for dom in domains.split(','):
        dom_nres = 0
        
        new_seg_chopping = []
        for segment in dom.split('_'):
            rng = segment.split('-')
            
            if len(rng) == 2:
                start, end = rng
                
                # Mainly for Chainsaw: 
                if args.offset_resi != 0:
                    start = str(int(start) + args.offset_resi)
                    end = str(int(end) + args.offset_resi)
                    
                seg_nres = int(end) - int(start) + 1

                if seg_nres >= args.min_fragment_size:
                    dom_nres += seg_nres
                    new_seg_chopping.append('-'.join([start, end]))
                    
        if dom_nres >= args.min_dom_size:
            new_dom_chopping.append('_'.join(natsorted(new_seg_chopping)))
            new_ndom += 1
            
    new_dom_chopping = ",".join(natsorted(new_dom_chopping))
    if len(new_dom_chopping) == 0:
        new_dom_chopping = 'na'
        
    return new_ndom, new_dom_chopping

def main():
    parser = argparse.ArgumentParser(description="Process input file and filter domain information.")
    parser.add_argument("input_file", help="Input file path")
    parser.add_argument("-o", "--output_file", help="Output file path")
    parser.add_argument("--changed_file", type=str, required=False, default=None, help="Output file containing changed targets")
    parser.add_argument("--min_dom_size", type=int, default=25, help="Minimum domain size (default: 25)")
    parser.add_argument("--min_fragment_size", type=int, default=5, help="Minimum fragment size (default: 5)")
    parser.add_argument("--offset_resi", type=int, default=0, help="Add this value to the chopping residue ids")
    args = parser.parse_args()
    
    if args.changed_file is None:
        args.changed_file = os.path.splitext(args.output_file)[0] + '.changed.txt'

    processed_data = {}
    highchanged, medchanged, lowchanged = 0, 0, 0
    changed = []

    with open(args.input_file, 'r') as f:
        for i, line in enumerate(f):
            line = line.rstrip('\n')

            target, md5, nres, nhigh, nmed, nlow, chophigh, chopmed, choplow = line.split()
            
            # print("          ", target, nhigh, nmed, nlow, chophigh, chopmed, choplow)
            
            if chophigh != 'na':
                new_nhigh, new_chophigh = filter_domains(chophigh, args)
                
                if new_chophigh != chophigh:
                    highchanged = 1
                    nhigh = new_nhigh
                    chophigh = new_chophigh
                    # print("     High:", target, new_nhigh, nmed, nlow, new_chophigh, chopmed, choplow)
                
            if chopmed != 'na':
                new_nmed, new_chopmed = filter_domains(chopmed, args)
                
                if new_chopmed != chopmed:
                    medchanged = 1
                    nmed = new_nmed
                    chopmed = new_chopmed
                    # print("      Med:", target, nhigh, new_nmed, nlow, chophigh, new_chopmed, choplow)
                
            if choplow != 'na':
                new_nlow, new_choplow = filter_domains(choplow, args)
                
                if new_choplow != choplow:
                    lowchanged = 1
                    nlow = new_nlow
                    choplow = new_choplow
                    # print("      Low:", target, nhigh, nmed, new_nlow, chophigh, chopmed, new_choplow)
                    
            if highchanged or medchanged or lowchanged:
                changed.append(target)
                    
            # with open(args.output_file, 'a+') as fn:
            #     fn.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            #         target, md5, nres, nhigh, nmed, nlow, chophigh, chopmed, choplow, proteome
            #     ))
                    
            processed_data[target] = {
                'md5': md5,
                'nres': nres,
                'nhigh': nhigh,
                'nmed': nmed,
                'nlow': nlow,
                'chophigh': chophigh,
                'chopmed': chopmed,
                'choplow': choplow,
            }
            
            # if i == 100:
            #     break
            
    # for target, data in processed_data.items():
    #     print(target, data['md5'], data['nres'], data['nhigh'], data['nmed'], data['nlow'],
    #                                                data['chophigh'], data['chopmed'], data['choplow'], data['proteome'])
                                                   
    with open(args.output_file, 'w') as output_file:
        for target, data in processed_data.items():
            line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(target, data['md5'], data['nres'], data['nhigh'], data['nmed'], data['nlow'],
                                                   data['chophigh'], data['chopmed'], data['choplow'])
            output_file.write(line + '\n')
            
    with open(args.changed_file, 'w') as fn:
        for target in changed: 
            fn.write(target + '\n')

if __name__ == "__main__":
    main()

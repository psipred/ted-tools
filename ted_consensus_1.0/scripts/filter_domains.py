import argparse

# This file is a part of TED: The Encyclopedia of Domains. If you utilize or reference any content from this file, 
# please cite the following paper:

# Lau et al., 2024. Exploring structural diversity across the protein universe with The Encyclopedia of Domains.

def main():
    parser = argparse.ArgumentParser(description="Process input file and filter domain information.")
    parser.add_argument("input_file", help="Input file path")
    parser.add_argument("-o", "--output_file", help="Output file path")
    parser.add_argument("--min_dom_size", type=int, default=25, help="Minimum domain size (default: 25)")
    parser.add_argument("--min_fragment_size", type=int, default=5, help="Minimum fragment size (default: 5)")
    parser.add_argument("--offset_resi", type=int, default=0, help="Add this value to the chopping residue ids")
    args = parser.parse_args()

    processed_data = {}  
    processed_domains = {}

    with open(args.input_file, 'r') as f:
        for i, line in enumerate(f):
            line = line.rstrip('\n')

            target, md5, nres, ndom, chopping, score = line.split()
            if target == 'chain_id':
                continue

            new_nres = 0
            new_chopping = []
            if chopping != 'NULL' and chopping != 'NO_SS':
                for nd, d in enumerate(chopping.split(',')):

                    dom_nres = 0
                    new_dd_chopping = []

                    for dd in d.split('_'):
                        rng = dd.split('-')
                        if len(rng) == 2:
                            start, end = rng
                            
                            # Mainly for Chainsaw: 
                            if args.offset_resi != 0:
                                start = str(int(start) + args.offset_resi)
                                end = str(int(end) + args.offset_resi)
                                
                            seg_nres = int(end) - int(start) + 1

                            if seg_nres >= args.min_fragment_size:
                                dom_nres += seg_nres
                                new_dd_chopping.append('-'.join([start, end]))

                    if dom_nres >= args.min_dom_size:
                        new_nres += dom_nres
                        
                        new_chopping.append('_'.join(new_dd_chopping))
                        
                        new_domain_id = f"{target}_{nd:02}"
                        
                        processed_domains[new_domain_id] = {
                            'md5': md5,
                            'nres': dom_nres,
                            'score': '1.000', # DUMMY
                            'plddt': '1.000', # DUMMY
                            'chopping': '_'.join(new_dd_chopping),
                        }

                new_ndom = len(new_chopping)
                new_chopping = ','.join(new_chopping)

                if len(new_chopping) == 0:
                    new_chopping = 'NULL'
            else:
                new_ndom = ndom
                new_chopping = chopping

            processed_data[target] = {
                'md5': md5,
                'nres': nres,
                'new_ndom': new_ndom,
                'new_chopping': new_chopping,
                'score': score
            }

    with open(args.output_file, 'w') as output_file:
        for target, data in processed_data.items():
            line = "{}\t{}\t{}\t{}\t{}\t{}".format(target, data['md5'], data['nres'], data['new_ndom'], data['new_chopping'], data['score'])
            output_file.write(line + '\n')

if __name__ == "__main__":
    main()

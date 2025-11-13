import os
import argparse
import subprocess
import hashlib
import re
import time

from natsort import natsorted

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
BINDIR = os.path.join(SCRIPT_DIR, 'bin')
UNIDOC = os.path.join(BINDIR, 'UniDoc_struct')
STRIDE = os.path.join(BINDIR, 'stride')
pdb_to_fasta = "pdb_tofasta"
pdb_selres= "pdb_selres"

resndict = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
            'PAD': 'X'}

def main():
    parser = argparse.ArgumentParser()
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--input", type=str, nargs="+", required=False, help="Pass a single file or list of filenames as -i 1ubq.pdb or *.pdb.")
    group.add_argument("-l", "--list", type=str, required=False, help="Pass a file containing paths to files.")
    
    parser.add_argument("-c",dest='chain', required=False, type=str, default='A', help="the chain of parsed protein")
    parser.add_argument("--out",dest='outfile', required=True, type=str, help="output file to write results to")
    parser.add_argument("--inherit_chopping", type=str, required=True, default=None, help="Pass a file containing choppings from Merizo or Chainsaw. Rows should match targets in --input.")
    
    args = parser.parse_args()
    
    if args.input is not None:
        files = args.input

    if args.list is not None:
        with open(args.list, 'r') as fn:
            files = [line.rstrip('\n') for line in fn]
            
    chopping_dict = {}
    with open(args.inherit_chopping, 'r') as fn:
        for line in fn:
            line_split = line.rstrip('\n').split('\t')
            chopping_dict[line_split[0]] = line_split[-2]

    with open(args.outfile, 'w') as fn:
        for pdb_path in files:
            start_time = time.time()
            bn, ext = os.path.splitext(os.path.basename(pdb_path))
            chopping = chopping_dict[bn]

            pdb = os.path.realpath(pdb_path)
            pdb_bn, pdb_ext = os.path.splitext(pdb)
            pdb_ss = pdb + '.ss'
            
            fasta = str(subprocess.check_output(f"{pdb_to_fasta} {pdb}", shell=True), 'utf-8').split('\n')[1:-1]
            seq = ''.join(fasta)
            
            md5 = hashlib.md5(seq.encode('utf-8')).hexdigest()
            nres = len(seq)
            
            # pdb_selres.py -1:10,20:30 pdb.pdb
            # If chopping is provided, then extract all domain residues from PDB using pdb_tools
            # save as new file and point to the new file
            if chopping not in ['NULL','NO_SS']:
                print(f"Processing chopping for {bn}: {chopping}")
                pdb_path_chopped = pdb_bn + '_chopped.pdb'
                resrng = chopping.replace('-',':').replace('_',',') # convert domain chopping into segments
                print(f"Extracting residues {resrng} from {pdb_path}")
                subprocess.check_output(f"{pdb_selres} -{resrng} {pdb_path} > {pdb_path_chopped} 2> /dev/null", shell=True)
                
                is_empty = os.stat(pdb_path_chopped).st_size == 0
                if_exist = os.path.exists(pdb_path_chopped)
                if if_exist and not is_empty:
                    print(f"Successfully created chopped PDB file: {pdb_path_chopped}")
                    pdb = pdb_path_chopped
                else:
                    print(f"Warning: Chopped PDB file is empty or does not exist: {pdb_path_chopped}")
                    
                try:
                    # Run secondary structure calculation with STRIDE
                    # print(f"Running STRIDE on {pdb} for chain {args.chain}")
                    # subprocess.check_output(f"{STRIDE} {pdb} -r{args.chain} > {pdb_ss} 2> /dev/null", shell=True)
                    
                    # Run UniDoc
                    # UniDoc_struct needs to be run from the directory containing it as it looks for ./stride
                    # Change to the directory containing UNIDOC
                    unidoc_dir = os.path.dirname(UNIDOC)
                    original_dir = os.getcwd()
                    os.chdir(unidoc_dir)
                    print(f"Running UniDoc on {pdb} in directory {unidoc_dir}")
                    
                    # Run UniDoc from its directory
                    output = subprocess.check_output(f"./UniDoc_struct {pdb} {args.chain}", shell=True)
                    
                    # Change back to original directory
                    os.chdir(original_dir)
                    # Format the output
                    output = str(output, 'utf-8').replace('~','-').replace(',','_').replace('/',',').rstrip('\n')
                    print(f"UniDoc output: {output}")

                    domains = output.split(',')
                    ndoms = len(domains)
                    chopping = ','.join(natsorted(domains))
                    print(f"Found {ndoms} domains: {chopping}")
                    
                    if chopping == '':
                        print("No domains found, setting chopping to NULL")
                        chopping = "NULL"
                        ndoms = 0
                    
                except Exception as e:
                    print(f"Error processing {pdb}: {str(e)}")
                    chopping = 'NO_SS'
                    ndoms = 0
            else:
                print(f"Skipping processing for {bn} - chopping is {chopping}")
                pdb_path_chopped = None
                ndoms = 0

            # end_time = time.time() - start_time 
            
            fn.write("{}\t{}\t{}\t{}\t{}\t{:.5f}\n".format(
                bn,
                md5,
                nres, 
                ndoms, 
                chopping,
                1.,
                # end_time,
            ))
            
            # Cleanup temporary files
            try:
                if os.path.exists(pdb_ss):
                    os.remove(pdb_ss)
                    print(f"Cleaned up temporary file: {pdb_ss}")

                if pdb_path_chopped and os.path.exists(pdb_path_chopped):
                    os.remove(pdb_path_chopped)
                    print(f"Cleaned up temporary file: {pdb_path_chopped}")
            except Exception as e:
                print(f"Warning: Error during cleanup: {str(e)}")
            
            print("--------------------------------------------------------------")
    
if __name__ == "__main__":
    main()
#!/usr/bin/env python

# Foldclass dbsearch (C) 2024 University College London

import os
import argparse
import pickle
import re
import subprocess

from math import sqrt, log

import numpy as np

import torch
import torch.nn as nn
import torch.nn.functional as F

from nndef_fold_egnn_embed import FoldClassNet


def run_tmalign(structure1_path, structure2_path, options=None):
    # Path to the tmalign executable (assume it is in current PATH)
    tmalign_path = 'tmalign'

    # Run tmalign as a subprocess
    if options is None:
        process = subprocess.Popen([tmalign_path, structure1_path, structure2_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    else:
        process = subprocess.Popen([tmalign_path, structure1_path, structure2_path, options], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    output, error = process.communicate()

    if process.returncode != 0:
        print(f"Error running tmalign: {error}")
        return None

    return output


def extract_tmalign_values(tmalign_output):
    # Define regular expressions to extract values
    aligned_length_pattern = re.compile(r'Aligned length=\s*(\d+),\s+RMSD=\s*([0-9.]+),\s+Seq_ID=n_identical/n_aligned=\s*([0-9.]+)')
    tm_score_pattern = re.compile(r'TM-score=\s*([0-9.]+)')
    
    # Extract values using regular expressions
    aligned_length_match = aligned_length_pattern.search(tmalign_output)
    tm_score_matches = tm_score_pattern.finditer(tmalign_output)

    # Extract values
    aligned_length = int(aligned_length_match.group(1)) if aligned_length_match else None
    rmsd = float(aligned_length_match.group(2)) if aligned_length_match else None
    seq_identity = float(aligned_length_match.group(3)) if aligned_length_match else None
    tm_scores = [float(match.group(1)) for match in tm_score_matches]

    # Capture three lines of alignment
    alignment_start_index = tmalign_output.find('(":" denotes residue pairs')
    alignment_lines = tmalign_output[alignment_start_index:].split('\n')[1:4]

    return aligned_length, rmsd, seq_identity, tm_scores, alignment_lines

def main():

    parser = argparse.ArgumentParser()
    # Add arguments
    parser.add_argument('-d', '--device', type=str, default='cpu', required=False)
    parser.add_argument('-j', '--skipk', type=int, default=0, required=False)
    parser.add_argument('-k', '--topk', type=int, default=1, required=False)
    parser.add_argument('-f', '--fastmode', type=int, default=1, required=False)
    parser.add_argument('-l', '--firstonly', type=int, default=1, required=False)
    parser.add_argument('-b', '--printbest', type=int, default=0, required=False)
    parser.add_argument('-n', '--dbname', type=str, required=True)
    parser.add_argument('-t', '--threads', type=int, default=-1, required=False)
    parser.add_argument('-s', '--mincos', type=float, default=0.5, required=False)
    parser.add_argument('-m', '--mintm', type=float, default=0.5, required=False)
    parser.add_argument('-r', '--maxrmsd', type=float, default=1000.0, required=False)
    parser.add_argument('-c', '--mincov', type=float, default=0.7, required=False)
    parser.add_argument('-p', '--pdb_path', type=str, default=None, required=False)
    parser.add_argument('pdbfile', type=str, nargs='+', help='PDB format file to process')
    # Parse the argument
    args = parser.parse_args()

    device = torch.device(args.device)
    topk = args.topk
    skipk = args.skipk
    dbname = args.dbname
    tmopts = '-fast' if args.fastmode else None 

    if args.threads > 0:
        torch.set_num_threads(args.threads)

    # Create neural network model
    network = FoldClassNet(128).to(device).eval().to(device)

    scriptdir = os.path.dirname(os.path.realpath(__file__))

    network.load_state_dict(torch.load(scriptdir + '/FINAL_foldclass_model.pt', map_location=lambda storage, loc: storage), strict=False)

    search_tensor = torch.load(dbname + '.pt').to(device)
    with open(dbname + '.targets', 'rb') as targfile:
        targets = pickle.load(targfile)

    assert len(targets) == search_tensor.size(0)

    target_lengths = torch.tensor([len(target[2]) for target in targets], dtype=torch.float, device=device)

    aa_conversion = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    for fname_q in args.pdbfile:
        with open(fname_q, 'r') as pdbfile:
            coords = []
            seq = []
            n = 0
            for line in pdbfile:
                if line[:4] == 'ATOM' and line[12:16] == ' CA ':
                    # Split the line
                    pdb_fields = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
                    coords.append(np.array([float(pdb_fields[6]), float(pdb_fields[7]), float(pdb_fields[8])]))
                    seq.append(aa_conversion.get(pdb_fields[3], 'X'))

        ca_coords_q = np.asarray(coords, dtype=np.float32)[:2000]
        seq_q = ''.join(seq[:2000])

        inputs = torch.from_numpy(ca_coords_q).unsqueeze(0).to(device)

        with torch.no_grad():
            qvec = network(inputs)

        mask = (target_lengths < len(seq_q) * args.mincov).float()
        scores = F.cosine_similarity(search_tensor, qvec, dim = -1) - mask * 10

        if topk > 0:
            this_topk = topk
            top_scores, top_indices = torch.topk(scores, skipk+topk, dim=0)
        else:
            this_topk = int((1-mask).sum())
            top_scores, top_indices = torch.topk(scores, skipk+this_topk, dim=0)

        bestscore = 0
        outflag = False
        for i in range(skipk, min(this_topk, scores.size(0))):
            fname_t, ca_coords_t, seq_t = targets[top_indices[i]]
            if args.pdb_path is not None:
                print(args.pdb_path, os.path.splitext(os.path.basename(fname_t))[0])
                fname_t = os.path.join(args.pdb_path, os.path.splitext(os.path.basename(fname_t))[0])
            if top_scores[i] >= args.mincos and len(seq_t) / len(seq_q) >= args.mincov:
                tmalign_output = run_tmalign(fname_q, fname_t, options=tmopts)
                if tmalign_output:
                    aligned_length, rmsd, seq_identity, tm_scores, alignment_lines = extract_tmalign_values(tmalign_output)
                    if aligned_length is not None and aligned_length >= len(seq_q) * args.mincov and max(tm_scores) >= args.mintm and rmsd < args.maxrmsd:
                        print(fname_q, fname_t, aligned_length, rmsd, seq_identity, len(seq_q), len(seq_t), '%.6f' % top_scores[i].item(), i, max(tm_scores))
                        outflag = True
                        if args.firstonly:
                            break
                    elif aligned_length is not None and aligned_length * max(tm_scores) > bestscore:
                        bestresult = (fname_q, fname_t, aligned_length, rmsd, seq_identity, len(seq_q), len(seq_t), '%.6f' % top_scores[i].item(), i, max(tm_scores))
                        bestscore = aligned_length * max(tm_scores)
        if args.printbest and not outflag and bestscore > 0:
            print(*bestresult)
if __name__=="__main__":
    main()

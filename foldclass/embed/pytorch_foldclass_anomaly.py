#!/usr/bin/env python

# By David T. Jones (C) 2024 University College London

import os
import argparse
import pickle

import numpy as np

import torch
import torch.nn as nn
import torch.nn.functional as F

from nndef_fold_egnn_embed import FoldClassNet


def main():

    parser = argparse.ArgumentParser()
    # Add arguments
    parser.add_argument('-d', '--device', type=str, default='cpu', required=False)
    parser.add_argument('-m', '--metric', type=int, default=0, required=False)
    parser.add_argument('-k', '--topk', type=int, default=50, required=False)
    parser.add_argument('-n', '--dbname', type=str, required=True)
    parser.add_argument('-t', '--threads', type=int, default=-1, required=False)
    parser.add_argument('pdbfile', type=str, nargs='+', help='PDB format file to process')
    # Parse the argument
    args = parser.parse_args()

    device = torch.device(args.device)
    topk = args.topk
    dbname = args.dbname

    if args.threads > 0:
        torch.set_num_threads(args.threads)

    # Create neural network model
    network = FoldClassNet(128).to(device).eval().to(device)

    scriptdir = os.path.dirname(os.path.realpath(__file__))

    network.load_state_dict(torch.load(scriptdir + '/FINAL_foldclass_model.pt', map_location=lambda storage, loc: storage), strict=False)

    search_tensor = torch.load(dbname + '.pt').to(device)

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

        if args.metric == 0:
            scores = (search_tensor - qvec).pow(2).sum(dim = -1).sqrt()
            top_scores, top_indices = torch.topk(scores, topk, largest=False, dim=0)
        else:
            scores = F.cosine_similarity(search_tensor, qvec, dim = -1)
            top_scores, top_indices = torch.topk(scores, topk, largest=True, dim=0)

        print(fname_q, top_scores.mean().item())

if __name__=="__main__":
    main()

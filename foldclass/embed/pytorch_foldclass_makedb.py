#!/usr/bin/env python

# By David T. Jones (C) 2024 University College London

import os
import argparse
import pickle
import tarfile

import numpy as np

import torch
import torch.nn as nn
import torch.nn.functional as F

from nndef_fold_egnn_embed import FoldClassNet

def main():

    parser = argparse.ArgumentParser()
    # Add arguments
    parser.add_argument('-d', '--device', type=str, default='cpu', required=False)
    parser.add_argument('-o', '--outfile', type=str, required=True)
    parser.add_argument('pdbfile', type=str, nargs='+', help='PDB format file to process')
    # Parse the argument
    args = parser.parse_args()

    device = torch.device(args.device)

    # Create neural network model
    network = FoldClassNet(128).to(device).eval().to(device)

    scriptdir = os.path.dirname(os.path.realpath(__file__))

    network.load_state_dict(torch.load(scriptdir + '/FINAL_foldclass_model.pt', map_location=lambda storage, loc: storage), strict=False)

    aa_conversion = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    targets = []
    tvecs = []
    coordset = set()

    def procfile(targpdbfile, fname):
        coords = []
        seq = []
        for line in targpdbfile:
            if isinstance(line, (bytes, bytearray)):
                line = line.decode('ascii')
            if line[:4] == 'ATOM' and line[12:16] == ' CA ':
                # Split the line
                pdb_fields = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
                coords.append(np.array([float(pdb_fields[6]), float(pdb_fields[7]), float(pdb_fields[8])]))
                seq.append(aa_conversion.get(pdb_fields[3], 'X'))

        if len(seq) < 5:
            return

        seq_t = ''.join(seq[:2000])

        ca_coords_t = np.array(coords, dtype=np.float32)[:2000]
        coordbytes = ca_coords_t.tobytes()

        if coordbytes in coordset:
            print("Skipped", fname)
            return
        else:
            coordset.add(coordbytes)

        inputs = torch.from_numpy(ca_coords_t).unsqueeze(0).to(device)

        with torch.no_grad():
            tvec = network(inputs).cpu()
            tvecs.append(tvec)
            #basename = os.path.splitext(os.path.basename(fname))[0]
            print(fname)
            targets.append((fname, ca_coords_t, seq_t))

    for fname in args.pdbfile:
        if tarfile.is_tarfile(fname):
            with tarfile.open(fname, 'r') as tar:
                for member in tar.getmembers():
                    if member.isfile():
                        procfile(tar.extractfile(member), member.name)
        else:
            procfile(open(fname, 'r'), fname)

    print("Number of target structures: ", len(tvecs))

    search_tensor = torch.cat(tvecs, dim=0)
    torch.save(search_tensor, args.outfile + '.pt')

    with open(args.outfile + '.targets', 'wb') as targfile:
        pickle.dump(targets, targfile)

if __name__=="__main__":
    main()

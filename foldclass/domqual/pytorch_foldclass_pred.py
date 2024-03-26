#!/usr/bin/env python

# By David T. Jones (C) 2024 University College London

from __future__ import print_function

import os
import argparse
import tarfile

import numpy as np

import torch
import torch.nn as nn
import torch.nn.functional as F

from nndef_fold_egnn import FoldClassNet


def main():

    parser = argparse.ArgumentParser()
    # Add arguments
    parser.add_argument('-d', '--device', type=str, default='cpu', required=False)
    parser.add_argument('pdbfile', type=str, nargs='+', help='PDB format file to process')
    # Parse the argument
    args = parser.parse_args()

    device = torch.device(args.device)

    # Create neural network model
    network = FoldClassNet(128, 3).to(device).eval().to(device)

    scriptdir = os.path.dirname(os.path.realpath(__file__))

    network.load_state_dict(torch.load(scriptdir + '/FINAL_foldclass_model.pt', map_location=lambda storage, loc: storage))


    def procfile(infile, fname):
        coords = []
        for line in infile:
            if line[:4] == 'ATOM' and line[12:16] == ' CA ':
                # Split the line
                pdb_fields = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
                coords.append(np.array([float(pdb_fields[6]), float(pdb_fields[7]), float(pdb_fields[8])]))

        length = len(coords)

        if length >= 5:
            ca_coords = np.asarray(coords, dtype=np.float32)[:length]

            inputs = torch.from_numpy(ca_coords).unsqueeze(0).to(device)

            pred = torch.softmax(network(inputs), dim=1)

            print(fname, pred[0,1].item())
        else:
            print(fname, 0.0)

    for fname in args.pdbfile:
        if fname.endswith('.tar.gz'):
            with tarfile.open(fname, 'r:gz') as tar:
                for member in tar.getmembers():
                    # Check if the file is a .pdb file
                    if member.isfile() and member.name.endswith('.pdb'):
                        procfile(tar.extractfile(member), member.name)
        else:
            procfile(open(fname, 'r'), fname)

if __name__=="__main__":
    main()

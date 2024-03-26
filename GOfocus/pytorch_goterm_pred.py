#!/usr/bin/env python

# GO term prediction with a dilated convolutional network and focal loss
# by David T. Jones (C) 2019 University College London

import sys
import os
import time
import random

from math import sqrt, log

from copy import deepcopy

import numpy as np

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim


#from nndef_goterm_transformer_cpu import TransformerNet
from nndef_goterm_pconv import PoolNet

# Read FASTA file
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(''.join(line.split()))
    if name: yield (name, ''.join(seq))


# ############################## Main program ################################
# Everything else will be handled in our main program now. We could pull out
# more functions to better separate the code, but it wouldn't make it any
# easier to read.

def main():
    golist = []

    setnum = sys.argv[2]

    NGOTERMS = 0
    with open('goterms_' + setnum + '.lst', 'r') as gotermfile:
        for line in gotermfile:
            goterm = line.rstrip().split()[0]
            golist.append(goterm)
            NGOTERMS += 1

    network = PoolNet(512,16,NGOTERMS).eval().cuda()

    network.load_state_dict(torch.load('set' + setnum + '/maxpool_goterm_model.pt', map_location=lambda storage, loc: storage))

    aa_trans = str.maketrans('ARNDCQEGHILKMFPSTWYVBJOUXZ', 'ABCDEFGHIJKLMNOPQRSTUUUUUU')

    with open(sys.argv[1]) as fp:
        for name, seq in read_fasta(fp):
            name = name[1:]
            seq = ''.join([i for i in seq if i.isalpha()]).upper()
            length = len(seq)
            seq = (np.frombuffer(seq.translate(aa_trans).encode('latin-1'), dtype=np.uint8) - ord('A')).reshape(1,length)

            inputs = torch.from_numpy(seq).type(torch.LongTensor).cuda()
            output = torch.sigmoid(network(inputs))

            for i, goterm in enumerate(golist):
                if output[0,i].item() > 0.63:
                    print("{}\t{}\t{:.3f}".format(name.split()[0], goterm, output[0,i].item()))

if __name__=="__main__":
    main()

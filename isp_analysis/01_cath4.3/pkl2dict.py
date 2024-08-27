#!/usr/bin/env python3

import pickle

def load_domdict(fname):
    p = pickle.load(open(fname, 'rb'))
    return p

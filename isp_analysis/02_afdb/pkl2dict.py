#!/usr/bin/env python3

import pickle
import numpy as np

def load_domdict(fname):
    p = pickle.load(open(fname, 'rb'))
    return p

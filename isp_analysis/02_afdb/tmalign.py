#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Adapted from https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/tmalign.py
(c) 2011-2012 Thomas Holder, MPI for Developmental Biology
LICENSE: BSD-2-Clause (https://opensource.org/license/bsd-2-clause/)

Modifications for PyMOL independence by Shaun Kandathil (SMK), Nov 2023

"""

import subprocess
# import tempfile
import os,sys
import re
import numpy as np

tmalign_binary = os.path.join(os.environ['HOME'], "software_staging", "tmalign", "tmalign")


def tmalign(mobile_filename, target_filename, args='', exe=tmalign_binary, tmpdir='/scratch0/skandath', quiet=True): #, ter=0, transform=1, object=None
    '''
DESCRIPTION

    TMalign wrapper returning the TM-score,
        resulting seq alignment, translation vector and rotation matrix

    Reference: Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9
    http://zhanglab.ccmb.med.umich.edu/TM-align/

ARGUMENTS

    mobile_filename, target_filename = string: paths to PDB files

    args = string: Extra arguments like -d0 5 -L 100 (DO NOT set -m or -outfmt)

    exe = string: Path to TMalign executable

    ter = 0/1: If ter=0, then ignore chain breaks because TMalign will stop
    at first TER record {default: 0}

    '''

    # SMK args `ter`, `transform`, and `object`are ignored in this version

    # ter, quiet = int(ter), int(quiet)

    # mobile_filename = tempfile.mktemp('.pdb', 'mobile')
    # target_filename = tempfile.mktemp('.pdb', 'target')

    # SMK This may not be a 'cluster-safe' option when using multiple processes
    # mfh, matrix_filename = tempfile.mkstemp('.txt', 'matrix')
    # os.close(mfh) # mkstemp() opens the fd, so close it

    mf, tf = [os.path.basename(fn) for fn in (mobile_filename, target_filename)]
    matrix_filename = os.path.join(tmpdir, '-'.join([mf, tf])+".mtx")

    errflag=0

    # mobile_ca_sele = '(%s) and (not hetatm) and name CA and alt +A' % (mobile)
    # target_ca_sele = '(%s) and (not hetatm) and name CA and alt +A' % (target)

    #exe = cmd.exp_path(exe)
    args = [exe, mobile_filename, target_filename, '-m', matrix_filename] + args.split()
    retry=3
    while retry>0:
        try:
            # process = subprocess.Popen(args, stdout=subprocess.PIPE, universal_newlines=True)
            process = subprocess.run(args, capture_output=True, universal_newlines=True)
            lines = process.stdout.split('\n')
            errlog = process.stderr.split('\n')
            if process.returncode != 0:
                errflag=1

            break
        except OSError as e:
            print(e, file=sys.stderr)
            print('Cannot execute "%s", retry...' % (exe), file=sys.stderr)
            retry-=1
    if retry <= 0:
        print('Fatal TMalign failure.', file=sys.stderr)
        sys.exit(1)

    # TMalign >= 2012/04/17
    if os.path.exists(matrix_filename):
        f = open(matrix_filename, 'r')
        lines += f.readlines()
        f.close()
        os.remove(matrix_filename)
    
    r = None
    re_score = re.compile(r'TM-score\s*=\s*(\d*\.\d*)')
    rowcount = 0
    # matrix = []
    t = []
    rot = []
    line_it = iter(lines)
    headercheck = False
    alignment = []
    for line in line_it:
        if 4 >= rowcount > 0:
            if rowcount >= 2:
                #a = list(map(np.float64, line.split()))
                a = line.strip().split()
                rot.extend(a[2:5])
                t.append(a[1])
            rowcount += 1
        elif not headercheck and line.startswith(' * '):
            a = line.split(None, 2)
            if len(a) == 3:
                headercheck = a[1]
        elif line.lower().startswith('------ the rotation matrix'):
            rowcount = 1
        elif line.startswith('(":" denotes'):
            alignment = [next(line_it).rstrip('\n') for i in range(3)]
        else:
            match = re_score.search(line)
            if match is not None:
                r = float(match.group(1))
        if not quiet:
            print(line.rstrip())

    if not quiet:
        for i in range(0, len(alignment[0]) - 1, 78):
            for line in alignment:
                print(line[i:i + 78])
            print('')

    #assert len(matrix) == 3 * 4
    try:
        assert len(t) == 3
        assert len(rot) == 9
    except AssertionError:
        print("ERROR aligning", mobile_filename, "and", target_filename, "- raw tm-align stderr:", file=sys.stderr)
        print(errlog, file=sys.stderr)
        errflag = 2
        rot = np.eye(3)
        t = np.zeros((3))
    #    matrix.extend([0.0, 0.0, 0.0, 1.0])

    # if int(transform):
    #     cmd.transform_selection('byobject (%s)' % (mobile), matrix, homogenous=1)

    # # alignment object
    # if object is not None:
    #     mobile_idx, target_idx = [], []
    #     space = {'mobile_idx': mobile_idx, 'target_idx': target_idx}
    #     cmd.iterate(mobile_ca_sele, 'mobile_idx.append("%s`%d" % (model, index))', space=space)
    #     cmd.iterate(target_ca_sele, 'target_idx.append("%s`%d" % (model, index))', space=space)
    #     for i, aa in enumerate(alignment[0]):
    #         if aa == '-':
    #             mobile_idx.insert(i, None)
    #     for i, aa in enumerate(alignment[2]):
    #         if aa == '-':
    #             target_idx.insert(i, None)
    #     if (len(mobile_idx) == len(target_idx) == len(alignment[2])):
    #         cmd.rms_cur(
    #             ' '.join(idx for (idx, m) in zip(mobile_idx, alignment[1]) if m in ':.'),
    #             ' '.join(idx for (idx, m) in zip(target_idx, alignment[1]) if m in ':.'),
    #             cycles=0, matchmaker=4, object=object)
    #     else:
    #         print('Could not load alignment object')

    if not quiet and r is not None:
        print('TM-align score = %.4f' % (r))

    return r, alignment, np.float64(t), np.float64(rot).reshape((3,3)), errflag


if __name__ == '__main__':

    from Bio.PDB import MMCIFParser, PDBParser, PDBIO

    basedir = os.path.dirname(__file__)

    model = os.path.join(basedir, "3nzqB02.pdb")
    native = os.path.join(basedir, "3sy1A00.pdb")

    tmscore, aln, transl, rot, errflag = tmalign(model, native, args="-ter 0 -o ~/TM_sup")

    p = PDBParser()
    m = p.get_structure('3nzqB02', model)
    n = p.get_structure('3sy1A00', native)

    for a in m.get_atoms():
        a.transform(np.transpose(rot), transl)

    io = PDBIO()
    io.set_structure(m)
    io.save(os.path.join(basedir, "3nzqB02_aligned.pdb"))

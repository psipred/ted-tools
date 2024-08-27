#!/usr/bin/env python3

# extract PDB files for 'master reference' domains used for alignment.

# Shaun Kandathil, 2023

import os, sys
import pickle

from subprocess import Popen, PIPE, run

def exec_(cmd):
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = proc.communicate()
    print(stderr)

refdom_d_fname = "/cluster/project9/afdb_domain_ext/results/dompair-eval-inputs/ref_domains_by_tedid.pkl"

refdom_d = pickle.load(open(refdom_d_fname, 'rb'))
# '3.10.450.140': {'taxid': '562-28', 'chain': 'AF-A0A376W4B5-F1-model_v4_TED02', 'plddt': '79.3871'}

# says tgz but really zip
tgz_dir = os.path.join(os.environ['HOME'], 'afdb_domain', 'afdb_domain_ext', 'results', 'proteome_consensus_domains_zip')
tgz_prefix = "proteome-tax_id-"
tgz_suffix = "_v4.consensus_domains" # .tar.gz added elsewhere
# proteome-tax_id-53326-6_v4.consensus_domains.tar.gz

out_dir = "/cluster/project9/afdb_domain_ext/results/dompair-eval-inputs/master-ref-domain-pdbs"
os.makedirs(out_dir, exist_ok=True)

# invert domids/domfiles and taxids
result = dict()

print('Collecting...')

for k in refdom_d.keys():
    taxid = refdom_d[k]['taxid']
    chain = refdom_d[k]['chain']

    # copy the entire tgz to tmpdir (scratch)
    tgz_base = tgz_prefix + taxid + tgz_suffix
    tgz_path     = tgz_dir + '/' + tgz_base + '.zip'
    # tgz_f_in_tmp = tmpdir  + "/" + tgz_base + '.tar.gz'
    # if not os.path.exists(tgz_f_in_tmp):
    #     shutil.copy(tgz_path, tgz_f_in_tmp)

    if taxid not in result.keys():
        # extract

        # result[taxid] = ['tar', '-C', out_dir, '-xf', tgz_path]
        result[taxid] = ['unzip', '-q', '-o', "-d", out_dir, tgz_path]  # unzip -q -o -d "${outdir}" "${zip_path}" "${domains[@]}"

    domfname = chain+'.pdb'
    result[taxid].append(domfname)


cmds = result.values()

for cmd in cmds:
    # print(cmd)
    process = run(cmd, capture_output=True, universal_newlines=True)
    if process.returncode != 0:
        print([a for a in process.stderr.split('\n')])
        sys.exit(1)

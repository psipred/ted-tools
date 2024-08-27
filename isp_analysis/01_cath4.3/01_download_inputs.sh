#!/bin/bash

domlist_url='ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/all-releases/v4_3_0/cath-classification-data/cath-domain-list-v4_3_0.txt'

wget $domlist_url

ln -sfv cath-domain-list-v4_3_0.txt cath-domain-list.txt

# we use a local copy of the CATH domain PDB files; if instead you need to download the domain PDBs, use the provided script (with 8 threads) like this:

# bash download_domain_pdbs_from_cath.sh -i cath-domain-list.txt -d $HOME/cath-dompdb -j 8

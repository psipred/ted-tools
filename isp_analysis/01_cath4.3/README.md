# Interacting superfamily pairs analysis for CATH domains


## Description of the scripts:
- `01_download_inputs.sh`: Download CATH 4.3 classification data
- Optionally get CATH domain PDBs; see below
- `02_make-cath-domain-list-pickle.py`: Make a lookup structure of which domains are in which chains.
- `03_dompair-eval-main.py`: Conduct the domain pair analysis.

## CATH domain PDB files:
We use a local copy of the CATH domain PDB files; if instead you need to download the domain PDBs, use the provided script (with 8 threads) like this, AFTER running step 01 above:

```
bash download_domain_pdbs_from_cath.sh -i cath-domain-list.txt -d $HOME/cath-dompdb -j 8
```


# AI-based Protein Domain Segmentation for TED: The Encyclopedia of Domains

## Summary

Version 1.0. 
Updated 26/08/2024.

This directory is part of TED: The Encyclopedia of Domains, and contains example scripts for performing domain segmentation on AFDB models to yield TED consensus domains. 

If you find this resource useful or use programs and scripts in TED, please cite the following paper: 
Lau et al., 2024. Exploring structural diversity across the protein universe with The Encyclopedia of Domains. 

Contents:
```
domain_consensus
├── README.md
├── run_segmentation.sh
|
├── example
│   └── AF-A0A1V6M2Y0-F1-model_v4.pdb
|
├── scripts
│   ├── domain_consensus.py
│   ├── filter_domains.py
│   ├── filter_domains_consensus.py
│   ├── get_consensus.py
│   ├── score_utils.py
│   └── segment.sh
|
└── programs/
    ├── chainsaw
    ├── Merizo
    └── UniDoc
```

The main program that runs domain segmentation and consensus calculation is `run_segmentation.sh`.
The `example` directory contains an example AFDB model that can be used for testing. 
The `example_output` shows the expected outputs from running `run_segmentation.sh`.
`scripts` contains scripts needed for segmentation and consensus calcualtion.
`programs` contains the source code for the three different segmentation programs used in TED.

## Installation

Dependencies can be installed by running `bash setup.sh`. 

## Running the program

TED consensus domains can be calculated by passing a directory containing PDB files to `run_segmentation.sh` as: 

```
bash run_segmentation.sh -i example -o example_output
```

This will generate the `example_output` directory which will contain the following outputs:
```
example_output
├── chopping_chainsaw.txt
├── chopping_merizo.txt
├── chopping_unidoc.txt
├── chopping_chainsaw.log
├── chopping_merizo.log
├── chopping_unidoc.log
├── consensus.tsv
├── consensus.log
└── consensus.tsv.changed.txt
```

Output files named `chopping_*.txt` contain the individual outputs from Chainsaw, Merizo and UniDoc. Note that as UniDoc does not classify non-domain residues (NDRs), the input to UniDoc includes the output of Merizo to remove NDRs prior to segmentation:
1. PDB -> Merizo -> Merizo output
2. PDB + Merizo output -> UniDoc -> UniDoc output
3. PDB + Chainsaw -> Chainsaw output

`consensus.tsv` contains the consensus domains assigned from the outputs of the three segmentation outputs. As per the TED publication, consensus domains are classified into high, medium and low consensus levels depending on whether there is agreement between three (high) or two (medium) methods. Agreement is determined as 70% overlap between the domain ranges. Any domains that do not agree between methods is sorted into the low consensus category and do not form part of TED. The `consensus.tsv.changed.txt` file will list any parses that have been changed as a result of the filtering. By default, after consensus calculation, domains with fewer than 25 residues and segments fewer than 5 residues are filtered out.

`.log` files contain any stdout and stderr output from each step. 

The output format of `consensus.tsv` is: 
```
AF-A0A1V6M2Y0-F1-model_v4	17c1e5cdd97ce14fb41a72ac2f0a15a8	1272	1	2	6	1179-1270	97-134_198-301_350-375_401-422,454-1126	102-135,202-315_352-420,423-558,559-752,753-945,946-1132
```

Fields correspond to:
1. AFDB target ID
2. MD5 string of the sequence - this is the sequence of the target converted into an MD5 string and can be used to find duplicate sequences across the AFDB
3. Number of residues in the target
4. Number of high consensus domains
5. Number of medium consensus domains
6. Number of low consensus domains
7. Domain string of high consensus domains
8. Domain string of medium consensus domains
9. Domain string of low consensus domains

Domains are separated by commas (`,`) and discontinuous domain ranges are separated by underscores (`_`).
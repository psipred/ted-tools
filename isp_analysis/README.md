# Interacting Superfamily Pair (ISP) analysis scripts for TED and CATH

Author: Shaun M Kandathil

Contact: s.kandathil@ucl.ac.uk

## Introduction

These directories contain scripts used to conduct the analyses of Interacting Superfamily Pairs (ISPs) in the paper. Each has its own README.

## Requirements:

- TM-align (we used the C++ version from [here](https://zhanggroup.org/TM-align/))
- Python >= 3.10
  - biopython >= 1.80
    - numpy < 2.0.0
    - tqdm
- R >= 4.0.0
  - beanplot
  - reticulate

## Running:

The scripts in the directories are meant to be run in their numbered order. Scripts without numbering in their filenames are called by the numbered scripts.

> [!IMPORTANT]
> You will need to modify paths specified in the scripts. In particular, you will need to specify the path to a local TM-align executable.

The final analyses comparing CATH and TED are carried out in the last steps of the 'afdb' pipeline.

The scripts assume that you have local copies for all input files. There are comments within each script that should help. For the CATH domain PDBs, there is a script provided which will download them.

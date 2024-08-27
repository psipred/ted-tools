# ISP analysis on TED domains

This is an archive of the code used to carry out the ISP analysis on the TED domains.

In order to run these scripts to completion, you will need the following:

- TED-100 domain PDB files.
  - We cannot readily make these available on account of their size. When running this analysis we created zip (not .tar.gz nor .gz, but .zip) files of TED domains. Each zip file contains the PDBs from a single proteome/taxid shard as defined in the AFDB. These files are named similar to `proteome-tax_id-53326-6_v4.consensus_domains.zip`.
  - PDB files can be reconstructed using the data on the Zenodo page linked in our paper. You will need full-chain AFDB PDB files in order to create the domain PDBs.
  
- PAE JSON files for all AFDB entries (optional; we provide a version of the relevant script that fetches them from the web).

## Description of the scripts:

- `01_make-ted-domain-proteome-dict.py`: Using a domain-level TSV file containing information about TED domains in multi-domain structures only, make a dict that indexes proteomes, chains and domains.
- `02_extract-master-ref-domains.py`: Define so-called master reference domains for alignment and save the info.
- `03_dompair-eval-main.py`: Main domain pair evaluation script. Meant to be run in parallel on a cluster.
- `04_get_all_ids_with_isps.py`: Use the outputs from 03 to get a list of AFDB chains having interacting domains.
- `05_get_isp_proteome_list_from_afdb_ids.py`: Use the outputs from 04 to make a few additional lookup data structures.
- `06_update_iv_chunks_with_choppings.py`: Update the outputs from 03 to include domain choppings, needed to evaluate PAE scores between them
- `07_dompair-eval-main.pt2.fromtars.py`: Update the outputs from 06 to add PAE scores. Meant to be run in parallel on a cluster. This version uses a local copy of AFDB .tar files. Alternatively, `extra/07_dompair-eval-main.pt2.fromweb.py` will retrieve PAE files from the web - use with caution and consult your cluster admin for advice.
- `08_collect-iv-outputs.py`: Aggregate the outputs from 07 into one, large pickled dictionary.
- `09_AnalyseResults.R`: R script that takes the output from 08 to conduct final analysis, including calculating CIO scores, and making various plots.


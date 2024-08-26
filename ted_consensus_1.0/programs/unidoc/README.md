# UniDoc

Brief introduction for domain parsing by UniDoc

## Structure-based domain parsing

***Step 1:*** Calculate the secondary sturcture by STRIDE(http://webclu.bio.wzw.tum.de/stride/), we provide in `./bin/stride`

​			Example of the  secondary sturcture  is  *./example/1a4iA_ss*

***Step2:*** Parse domain using UniDoc, we provide in `./bin/UniDoc_structure`

​			**Usage:**`./bin/UniDoc_structure pdb chain pdb_ss`

​			Example: `./bin/UniDoc_structure ./example/1a4i.pdb A ./example/1a4iA_ss`

**Note:** You can use *./Run_UniDoc_from_scratch_structure.py* to  parse the domain from scratch

​			**Usage:**  `python Run_UniDoc_from_scratch_structure.py -i seq.pdb -c Chain  `

​			**mandatory arguments:**
				-i native structure in PDB format.
				-c the chain of native structure

​			**optional arguments:**

​					-o the output file

## Sequence-based domain parsing

***Step 1:***   Using the trRosettaX to obtain the inter-residue distance from the trRosettaX web server(https://yanglab.nankai.edu.cn/trRosetta/)  or the  the new version of trRosettaX software package. (Only the predicted distance is required) . 

​			Example of the inter-residue distance  is  *./example/1a04a.npz*

​			When the user get the prediction of the inter-residue distance, they should use the *./bin/npz2dist.py* to get the inter-residue distance matrix.

​			**Usage:** `python npz2dist.py -n seq.npz -o seq.dist.mat  `

​			Example of the inter-residue distance  is  *./example/1a04a.dist.mat*

***Step 2:***  Using PSIPRED to predict the secondary structure of target. (http://bioinf.cs.ucl.ac.uk/psipred/) or download the latest version of PSIPRED software package from http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/.  When the user get the prediction of the secondary structure, they should use the *./bin/format_ss.pl* to transform the secondary structure format .

​			**Usage:** `./bin/format_ss.pl seq.ss2 seq.horiz` (seq.ss2 and seq.horiz are from PSIPRED. More about them you can see the http://bioinf.cs.ucl.ac.uk/psipred/ )

​			Example of the secondary structure format  is  *./example/seq.dat*

***Step 3***: Parse domain using UniDoc, we provide in `./bin/UniDoc_sequence`

​			**Usage:**`./bin/UniDoc_sequence seq.dist.mat seq.dat`

​			Example:` ./bin/UniDoc_sequence ./example/1a04a.dist.mat ./example/seq.dat`



**Note:** You can use *./Run_UniDoc_from_scratch_sequence.py*  to  parse the domain from scratch, if you prepare  the the inter-residue distance(*seq.npz*) and the secondary structure from PSIPRED( *seq.ss2 and seq.horiz*)

​			**Usage:**  `python Run_UniDoc_from_scratch_sequence.py -f seq.ss2 -s seq.horiz -n seq.npz  `

​			**mandatory arguments:**

						-f the *.ss2 file from PSIPRED.
​						-s the *.horiz file from PSIPRED.
​						-n the npz file

​			**optional arguments:**

​					-o the output file





​	






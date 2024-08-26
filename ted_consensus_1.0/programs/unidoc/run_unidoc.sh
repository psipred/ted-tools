#!/usr/bin/env bash

# Author:
#   A Lau - 2023-05-24

# This script runs the UniDoc method on a provided list of PDB files.

# Usage: run_unidoc.sh <target_list.txt> <model_dir> <results.txt>

# If results.txt exists, will ask to overwrite. 

set -eu

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
unidoc="${SCRIPT_DIR}/Run_UniDoc_from_scratch_structure_afdb.py"
pdb_tofasta=$(which pdb_tofasta)

input=$(readlink -f $1)
model_dir=$(readlink -f $2)
output=$(readlink -f $3)

dom_delim=","   # character to count for ndom

if [ -e "$output" ]; then
    read -p "Output file '${output}' already exists. Overwrite? (y/n) " choice
    case "$choice" in
        y|Y )
            rm $output
            ;;
        n|N )
            exit 0
            ;;
        * )
            echo "Invalid choice. Exiting."
            exit 1
            ;;
    esac
fi

while read -r line; do
    pdb_path="${model_dir}/${line}.pdb"
    seq=$($pdb_tofasta $pdb_path | tail -n +2 | tr -d '\n')

    if test -f $pdb_path; then
        out=$(python $unidoc -i $pdb_path -c A -s $seq)
        chopping=$(cut -d " " -f1 <<< $out)
        nres=$(cut -d " " -f2 <<< $out)
        md5=$(cut -d " " -f3 <<< $out)

        ndelim="${chopping//[^$dom_delim]}"
        ndom="$(( ${#ndelim} +1 ))"
        echo -e "${line}\t${md5}\t${nres}\t${ndom}\t${chopping}\t1" >> ${output} # tab delimited

    else
        echo "${line}\t${md5}\tSegmentation Failed." # tab delimited
    fi

done < ${input}

echo "Done." >> ${output}

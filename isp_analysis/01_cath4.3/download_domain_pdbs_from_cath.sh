#!/bin/bash

# Andy's script to download domain PDBs from CATH via REST API

input_file="cath-domain-list.txt"

while getopts "i:d:j:" option; do
    case "${option}" in
        i) input_file=${OPTARG};;
        d) output_dir=${OPTARG};;
        j) num_jobs=${OPTARG};;
    esac
done

if [ -z "$input_file" ] || [ -z "$output_dir" ] || [ -z "$num_jobs" ]; then
    echo "Usage: $0 -i <input_file> -d <output_dir> -j <num_jobs>"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p $output_dir

# Function to download PDB file
download_pdb() {
    domain_id=$1
    output_dir=$2
    wget -q "http://www.cathdb.info/version/v4_3_0/api/rest/id/$domain_id.pdb" -O "$output_dir/$domain_id"
    echo "Downloaded $domain_id"
}

# Export the function to make it accessible by parallel
export -f download_pdb

# Read input file and start parallel downloads
grep -v '^#' $input_file | cut -d ' ' -f 1 | xargs -I {} -P $num_jobs bash -c 'download_pdb "$@"' _ {} $output_dir

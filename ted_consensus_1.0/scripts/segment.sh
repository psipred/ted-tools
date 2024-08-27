#!/bin/bash

# This file is a part of TED: The Encyclopedia of Domains. If you utilize or reference any content from this file, 
# please cite the following paper:

# Lau et al., 2024. Exploring structural diversity across the protein universe with The Encyclopedia of Domains.

# Script for running segmentation methods given a directory of structures.

# Usage: 
# bash run_segment_afdb.sh -i <structure_directory> -m <merizo/unidoc/chainsaw>

set -eu

# Directories and paths
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PROG_DIR="${SCRIPT_DIR}/../programs"

py=$(which python)
custom_chopping=''

FILTER_DOMAINS="${SCRIPT_DIR}/filter_domains.py"

while getopts ":i:m:o:c:" opt; do
  case $opt in
    i) inputs=$(readlink -f "$OPTARG") ;;
    m) method=$OPTARG ;;
    o) output=$(readlink -f "$OPTARG") ;;
    c) custom_chopping=${OPTARG} ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Check if both options are provided
if [[ -z "${inputs}" || -z "${method}" || -z "${output}" ]]; then
  echo "Usage: run_segment_afdb.sh -i <structure_directory> -m <merizo/unidoc/chainsaw> -o <output_directory> [-c <chopping>]"
  exit 1
fi

case $method in
  "merizo")
    RUN_SCRIPT="${PROG_DIR}/merizo/predict_afdb.py"
    OFFSET_RESI=0
    ;;
  "unidoc")
    RUN_SCRIPT="${PROG_DIR}/unidoc/Run_UniDoc_from_scratch_structure_afdb.py"
    OFFSET_RESI=0
    ;;
  "chainsaw")
    RUN_SCRIPT="${PROG_DIR}/chainsaw/get_predictions.py"
    OFFSET_RESI=1
    ;;
  *)
    echo "Invalid method: ${method}. Allowed options are 'merizo', 'unidoc', or 'chainsaw'."
    exit 1
    ;;
esac

# Run method
output_file="${output%/}/chopping_${method}.txt"

echo "Running ${method} on targets in ${inputs}"

# Each method will take the list containing the paths to the targets
if [ "${method}" = "merizo" ] || [ "${method}" = "unidoc" ]; then
    target_list="${output%/}targets.txt"
    readlink -f "${inputs}/"*.pdb > "${target_list}"

    if [[ ${custom_chopping} == '' ]]; then
        ${py} "${RUN_SCRIPT}" -l "${target_list}" --out "${output_file}"
    else
        ${py} "${RUN_SCRIPT}" -l "${target_list}" --out "${output_file}" --inherit_chopping "${custom_chopping}"
    fi

    # Cleanup
    if test -f "${target_list}"; then
        rm "${target_list}"
    fi

elif [ "${method}" = "chainsaw" ]; then
    ${py} "${RUN_SCRIPT}" --structure_directory "${inputs}" -o "${output_file}" --append
fi

# Filter choppings to remove small segments and single-residue domains
if test -f "${output_file}"; then
    "${py}" "${FILTER_DOMAINS}" "${output_file}" -o "${output_file}.tmp" --offset_resi "${OFFSET_RESI}"

    if [ $? == 0 ]; then
        mv "${output_file}.tmp" "${output_file}"
    fi
else
    echo "Expected to find output file at ${output_file}"
    exit 1
fi

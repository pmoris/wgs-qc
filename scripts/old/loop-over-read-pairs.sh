#!/usr/bin/env bash

##########################################################
# overview of best practices for looping over read pairs #
##########################################################

# set bash strict mode
set -euo pipefail
# equivalent to:
# set -o errexit  # exit script when a command fails
# set -o nounset  # exit script when unset variable is accessed
# set -o pipefail # fail pipeline when any of its commands fail

# allow debug mode by running `TRACE=1 ./script.sh` - equivalent to `set -x`
if [[ "${TRACE-0}" == "1" ]]; then set -o xtrace; fi

# get file path of script
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

# set number of threads for downstream tools
n_threads=8

# define in and outputs relative to script
fastq_dir="${SCRIPT_DIR}/../data/fastq/pf"
ref="${SCRIPT_DIR}/../data/ref/Pf3D7_01_v3.fa"
read_suffix="_R1_001.fastq.gz"
output_dir="${SCRIPT_DIR}/../output/bwa"
# mkdir -p "${output_dir}"

# check if fastq directory exist
if [ ! -d "${fastq_dir}" ]; then
    echo "FASTQ input directory (${fastq_dir}) does not exist."
fi

# check if reference fasta file exists
if ! [ -f "${ref}" ]; then
    echo "Reference fasta file not found (${ref})."
fi

# log run options
printf "
BWA alignment script | $(basename "$0")
==============================================

FASTQ reads directory:      ${fastq_dir}
Reference genome:           ${ref}
threads:                    ${n_threads}
Output directory:           ${output_dir}
"

###################
# start of script #
###################

# create reference index if it does not yet exist
for i in "${ref}."{amb,ann,bwt,pac,sa}; do
    if ! [ -f "${i}" ]; then
        bwa index "${ref}"
        break
    fi
done

# loop through fastq read pairs

# basename option
for read_1 in "${fastq_dir}/"*"${read_suffix}"; do
    sample_name=$(basename "${read_1}" "${read_suffix}")
    # alternatively, first strip suffix using ${read_1%_R1_001.fastq.gz}
    # and then use basename without the second argument
    read_2="${fastq_dir}/${sample_name}_R2_001.fastq.gz"
    echo sample="$sample_name" R1="$read_1" R2="$read_2" output="${output_dir}/${sample_name}.out"
done

## parameter expansion option
for read_1 in "${fastq_dir}/"*"${read_suffix}"; do
    # remove suffix to retrieve file path of sample to which suffix can be added
    sample_path="${read_1%"${read_suffix}"}"
    # remove filepath up to filename to retrieve sample name that can be used to create new named output files
    sample_name="${sample_path##*/}"
    echo bwa R1="${sample_path}_R1_001.fastq.gz" R2="${sample_path}_R2_001.fastq.gz" OUT="${output_dir}/${sample_name}.out"

    # replace suffix by stripping pattern and adding manually
    read_2="${read_1%"${read_suffix}"}_R2_001.fastq.gz"
    echo R2 = "${read_2}"

    # replace suffix by replacing using parameter expansion
    read_2="${read_1/${read_suffix}/_R2_001.fastq.gz}"
    echo R2 = "${read_2}"

    # shorter version as above
    read_2="${read_1/R1/R2}"
    echo R2 = "${read_2}"
done

# for r1 in "${fastq_dir}"/pf*_R1_001.fastq.gz; do

#     # get filepath containing basename of each read pair
#     sample_path=${r1%_R1_001.fastq.gz}

#     # convert to basename of each read without the filepath prefix
#     sample_name="$(basename ${sample_path})"

#     # create output filepath with a directory for each read pair
#     output_prefix="${output_dir}/${sample_name}"

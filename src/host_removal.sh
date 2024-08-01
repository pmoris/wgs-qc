#! /bin/env bash

#################################################################
# Script to perform quality control and trimming of fastq reads #
#################################################################

# set bash strict mode
set -euo pipefail

# allow debug mode by running `TRACE=1 ./qc.sh`
if [[ "${TRACE-0}" == "1" ]]; then set -x; fi

# get file path of script
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
# Alternatively, use absolute paths, but this makes the script less portable.

#####################
# Options and paths #
#####################

# set number of threads for downstream tools
n_threads=8

# define in and outputs
# by defining everything here, there is no need to `cd` to directories first
fastq_dir="${SCRIPT_DIR}/../data/fastq/pf"
output_dir="${SCRIPT_DIR}/../results/"
fastq_screen_conf="${SCRIPT_DIR}/fastq-screen.conf"
ref_human="${SCRIPT_DIR}/../data/ref/GRCh38.chr21.fa.gz"
ref_pf="${SCRIPT_DIR}/../data/ref/Pf3D7_01_v3.fa.gz"

# create output directories
mkdir -p "${output_dir}/fastqc" \
    "${output_dir}/fastq-screen" \
    "${output_dir}/fastp" \
    "${output_dir}/multiqc"

# define fastq read suffix
read_1_suffix="_R1_001.fastq.gz"
read_2_suffix="_R2_001.fastq.gz"

# check if fastq directory exist
if [ ! -d "${fastq_dir}" ]; then
    echo "FASTQ input directory (${fastq_dir}) does not exist."
fi

# check if reference fasta file exists
for ref in ${ref_human} ${ref_pf}; do
    if ! [ -f "${ref}" ]; then
        echo "Reference fasta file not found (${ref})."
    fi
done

# log run options
printf "
FastQ QC script | $(basename "$0")
==============================================

Output directory:           ${output_dir}
FASTQ reads directory:      ${fastq_dir}
Reference human:            ${ref_human}
Reference Pfalciparum:      ${ref_pf}
threads:                    ${n_threads}
"

###############
# Start of QC #
###############

#!/usr/bin/env bash

# set bash strict mode
set -euo pipefail

# allow debug mode by running `TRACE=1 ./script.sh` - equivalent to `set -x`
if [[ "${TRACE-0}" == "1" ]]; then set -o xtrace; fi

# get file path of script
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

# set number of threads for downstream tools
n_threads=8

# define in and outputs
ref="${SCRIPT_DIR}/../data/ref/Pf3D7_01_v3.fa"
# ref_human="${SCRIPT_DIR}/../data/ref/GRCh38.chr21.fa.gz"
# ref_concat="${SCRIPT_DIR}/../data/ref/concat.fa.gz"
fastq_dir="${SCRIPT_DIR}/../data/fastq/pf"
# reads_human="${SCRIPT_DIR}/../data/fastq/human"
output_dir="${SCRIPT_DIR}/../output/bwa/pf"
mkdir -p "${output_dir}"

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
BWA MEM script | $(basename "$0")
==============================================

Output directory:           ${output_dir}
FASTQ reads directory:      ${fastq_dir}
Reference pfalciparum:      ${ref}
threads:                    ${n_threads}
"

###################
# start of script #
###################

# create reference index if it does not yet exist
# for ref in ${reference_pk_H} ${reference_pk_H} ${reference_human}; do
for i in "${ref}."{amb,ann,bwt,pac,sa}; do
    if ! [ -f "${i}" ]; then
        echo "Building BWA index for ${ref}..."
        bwa index "${ref}"
        break
    fi
done
# done

# map fastq read pairs
for r1 in "${fastq_dir}"/*_R1_001.fastq.gz; do

    # get filepath containing basename of each read pair
    sample_path="${r1%_R1_001.fastq.gz}"

    # convert to basename of each read without the filepath prefix
    sample_name="${sample_path##*/}"

    # create output filepath for each read pair
    output_prefix="${output_dir}/${sample_name}"
    # mkdir -p "${output_prefix}"

    printf "\nProcessing sample %s ...\n\n" "${sample_name}"

    # map to reference genome
    bwa mem \
        -t "${n_threads}" \
        -Y -K 100000000 \
        -R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA" \
        "${ref}" \
        "${sample_path}_R1_001.fastq.gz" \
        "${sample_path}_R2_001.fastq.gz" |
        # sort and compress to bam
        samtools sort -@ "${n_threads}" \
            -o "${output_prefix}.sort.bam"
done

#!/usr/bin/env bash

# TODO: set suffix as variable
# TODO: set threads to n_threads
# TODO: compare with existing scripts
# TODO: change basename into parameter expansion?
# TODO: split different tasks over different loops
# TODO: provide references as list that can be re-used by log run options and loops?
# TODO: clean up read file extension clean up somehow
# TODO: BWA MEM script | $(basename "${BASH_SOURCE[0]}") vs BWA MEM script | $(basename "$0")

#################################################################
# Script to perform quality control and trimming of fastq reads #
#################################################################

# set bash strict mode
set -euo pipefail

# allow debug mode by running `TRACE=1 ./qc.sh`
if [[ "${TRACE-0}" == "1" ]]; then set -x; fi

# get file path of project root to allow it to be run from any working directory
# Otherwise, you would need to make sure to call the script from within the directory where it is stored.
# Alternatively, use absolute paths, but this makes the script less portable.
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
PROJECT_ROOT=$(realpath "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)/../../")
echo "Project root = ${PROJECT_ROOT}"

#####################
# Options and paths #
#####################

# set number of threads for downstream tools
n_threads="${SLURM_CPUS_PER_TASK:-8}"

# define in and outputs
# by defining everything here, there is no need to `cd` to directories first
fastq_dir="${PROJECT_ROOT}/data/fastq/pf"
output_dir="${PROJECT_ROOT}/results/"

# config files
fastq_screen_conf="${PROJECT_ROOT}/config/fastq-screen.conf"
multiqc_conf="${PROJECT_ROOT}/config/multiqc_config.yaml"

ref_human="${PROJECT_ROOT}/data/ref/GRCh38.chr21.fa.gz"
ref_pf="${PROJECT_ROOT}/data/ref/Pf3D7_01_v3.fa.gz"

# define fastq read suffix
read_1_suffix="_R1_001.fastq.gz"
read_2_suffix="_R2_001.fastq.gz"
read_file_extension=".fastq.gz"

# create output directories
mkdir -p "${output_dir}/fastqc" \
    "${output_dir}/fastq-screen" \
    "${output_dir}/fastp" \
    "${output_dir}/multiqc"

# check if fastq directory exist
if [ ! -d "${fastq_dir}" ]; then
    printf "\nFASTQ input directory (${fastq_dir}) does not exist.\n"
    exit 1
fi

# check if reference fasta file exists
for ref in ${ref_human} ${ref_pf}; do
    if ! [ -f "${ref}" ]; then
        printf "\nReference fasta file not found (${ref}).\n"
        exit 1
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
# References:      ${ref_pf} \t ${ref_pv} \t ${ref_poc} \t ${ref_pow} \t ${ref_pm} \t ${ref_phix}

###############
# Start of QC #
###############

# create reference index if it does not yet exist
# required for fastq-screen
for ref in ${ref_human} ${ref_pf}; do
    for i in "${ref}."{amb,ann,bwt,pac,sa}; do
        if ! [ -f "${i}" ]; then
            index_files_found=0
            printf "\nBuilding BWA index for ${ref}...\n"
            bwa index "${ref}"
            break
        else
            index_files_found=1
        fi
    done
    if [ "$index_files_found" -eq 1 ]; then
        printf "\nFound BWA index files for ${ref}, skipping indexing step...\n"
    fi
done

# run initial quality control
# note that --threads option only works when providing
# multiple files, i.e. we need to use a glob instead of a loop
# (alternatively use parallel:
# `find *.fq | parallel -j 10 "fastqc {} --outdir ...` or find exec )
printf "\nRunning FastQC prior to trimming...\n"
fastqc \
    --threads "${n_threads}" \
    --outdir "${output_dir}/fastqc" \
    "${fastq_dir}/"*.fastq.gz

# run fastq-screen (threads option inherited by bwa/bowtie)
printf "\nRunning FastQ Screen...\n"
for read in "${fastq_dir}/"*.fastq.gz; do
    fastq_screen \
        --threads "${n_threads}" \
        --aligner bwa \
        --conf "${fastq_screen_conf}" \
        --outdir "${output_dir}/fastq-screen" \
        "${read}"
done

# trim reads using fastp (quality and adapters) using built-in multi-threading
# fastp can use up to 16 threads, but efficiency is higher around 2-4 (https://hpc.nih.gov/training/gatk_tutorial/preproc.html#preproc-trim) => using gnu parallel would be more efficient for many samples
fastp_threads=$(( ${n_threads} > 16 ? 16 : ${n_threads} ))

for read_1 in "${fastq_dir}/"*"${read_1_suffix}"; do
    sample_name=$(basename "${read_1}" "${read_1_suffix}")
    read_2="${fastq_dir}/${sample_name}${read_2_suffix}"
    printf "\nRunning fastp for ${sample_name}...\n"

    fastp \
        --in1 "${read_1}" \
        --in2 "${read_2}" \
        --out1 "${output_dir}/fastp/${sample_name}${read_1_suffix%${read_file_extension}}.trim.fastq.gz" \
        --out2 "${output_dir}/fastp/${sample_name}${read_2_suffix%${read_file_extension}}.trim.fastq.gz" \
        --json "${output_dir}/fastp/${sample_name}.trim.json" \
        --html "${output_dir}/fastp/${sample_name}.trim.html" \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --thread "${fastp_threads}"

done

# Additional options to consider:
        # --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        # --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
# --adapter_fasta                  specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file (string [=])
# --unpaired1                      for PE input, if read1 passed QC but read2 not, it will be written to unpaired1. Default is to discard it. (string [=])
# --unpaired2                      for PE input, if read2 passed QC but read1 not, it will be written to unpaired2. If --unpaired2 is same as --unpaired1 (default mode), both unpaired reads will be written to this same file. (string [=])
# --failed_out                     specify the file to store reads that cannot pass the filters. (string [=])
# --trim_poly_x                    enable polyX trimming in 3' ends.
# --merge                          for paired-end input, merge each pair of reads into a single read if they are overlapped. The merged reads will be written to the file given by --merged_out, the unmerged reads will be written to the files specified by --out1 and --out2. The merging mode is disabled by default.
# --cut_front, --cut_tail, --cut_right => might interfere with deduplication
# -q, --qualified_quality_phred      the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
# -l, --length_required              reads shorter than length_required will be discarded, default is 15. (int [=15])

# re-run qc after trimming
printf "\nRe-running FastQC after trimming...\n"
fastqc \
    --threads "${n_threads}" \
    --outdir "${output_dir}/fastqc/" \
    "${output_dir}/fastp/"*".trim.fastq.gz"

# aggregate results with multiQC
multiqc --force "${output_dir}" --config "${multiqc_conf}" --outdir "${output_dir}/multiqc"

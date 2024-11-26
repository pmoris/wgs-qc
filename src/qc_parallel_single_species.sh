#!/usr/bin/env bash

# TODO: set suffix as variable
# TODO: set threads to n_threads
# TODO: compare with existing scripts
# TODO: change basename into parameter expansion?
# TODO: split different tasks over different loops
# TODO: provide references as list that can be re-used by log run options and loops?
# TODO: clean up read file extension clean up somehow

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
PROJECT_ROOT=$(realpath "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)/../")
echo "Project root = ${PROJECT_ROOT}"

#####################
# Options and paths #
#####################

# read species from command line
species=$1

# set number of threads for downstream tools
n_threads="${SLURM_CPUS_PER_TASK:-8}"

# define in and outputs
# by defining everything here, there is no need to `cd` to directories first
# ! trailing slash is needed for `find`
fastq_dir="${PROJECT_ROOT}/data/fastq/"
output_dir="${PROJECT_ROOT}/results/"

# reference files
ref_human="${PROJECT_ROOT}/data/ref/human/Homo_sapiens.GRCh38.p14.GENCODE.release45/GRCh38.primary_assembly.genome.fa.gz"
ref_pf="${PROJECT_ROOT}/data/ref/Pfalciparum/PlasmoDB-release-68/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"
ref_pv="${PROJECT_ROOT}/data/ref/Pvivax/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.fasta"
ref_pm="${PROJECT_ROOT}/data/ref/Pmalariae/PlasmoDB-release-68/PlasmoDB-68_PmalariaeUG01_Genome.fasta"
ref_poc="${PROJECT_ROOT}/data/ref/Povale/PlasmoDB-release-68/PlasmoDB-68_PovalecurtisiGH01_Genome.fasta"
ref_pow="${PROJECT_ROOT}/data/ref/Povale/PlasmoDB-release-68/PlasmoDB-68_PovalewallikeriPowCR01_Genome.fasta"
ref_pk="${PROJECT_ROOT}/data/ref/Pknowlesi/PlasmoDB-release-68/PlasmoDB-68_PknowlesiH_Genome.fasta"
ref_phix="${PROJECT_ROOT}/data/ref/PhiX/PhiX-NC_001422.1.fasta"

if [[ "${species}" == "pf" ]]; then
    ref_plasmodium="${ref_pf}"
elif [[ "${species}" == "pv" ]]; then
    ref_plasmodium="${ref_pv}"
elif [[ "${species}" == "pm" ]]; then
    ref_plasmodium="${ref_pm}"
elif [[ "${species}" == "poc" ]]; then
    ref_plasmodium="${ref_poc}"
elif [[ "${species}" == "pow" ]]; then
    ref_plasmodium="${ref_pow}"
elif [[ "${species}" == "pk" ]]; then
    ref_plasmodium="${ref_pk}"
fi

# config files
fastq_screen_conf="${PROJECT_ROOT}/config/fastq-screen-${species}.conf"
multiqc_conf="${PROJECT_ROOT}/config/multiqc_config.yaml"

# check if fastq directory exist
if [ ! -d "${fastq_dir}" ]; then
    echo "FASTQ input directory (${fastq_dir}) does not exist."
    exit 1
fi

# check if reference fasta files exists
for ref in ${ref_human} ${ref_plasmodium} ${ref_phix}; do
    if ! [ -f "${ref}" ]; then
        echo "Reference fasta file not found (${ref})."
        exit 1
    fi
done

# create output directories
mkdir -p "${output_dir}/fastqc" \
    "${output_dir}/fastq-screen" \
    "${output_dir}/fastp" \
    "${output_dir}/multiqc"

# log run options
printf "
FastQ QC script | $(basename "$0")
==============================================

Output directory:           ${output_dir}
FASTQ reads directory:      ${fastq_dir}
References:                 ${ref_human} \t ${ref_plasmodium} \t ${ref_phix}
threads:                    ${n_threads}
"

###############
# Start of QC #
###############

# create reference index if it does not yet exist
# required for fastq-screen
for ref in ${ref_human} ${ref_plasmodium} ${ref_phix}; do
    for i in "${ref}."{amb,ann,bwt,pac,sa}; do
        if ! [ -f "${i}" ]; then
            index_files_found=0
            echo "Building BWA index for ${ref}..."
            bwa index "${ref}"
            break
        else
            index_files_found=1
        fi
    done
    if [ "$index_files_found" -eq 1 ]; then
        echo "Found BWA index files for ${ref}, skipping indexing step..."
    fi
done

# run initial quality control
# note that --threads option only works when providing
# multiple files, i.e. we need to use a glob instead of a loop
# (alternatively use parallel:
# `find *.fq | parallel -j 10 "fastqc {} --outdir ...` or find exec )
echo "Running FastQC prior to trimming..."
fastqc \
    --threads "${n_threads}" \
    --outdir "${output_dir}/fastqc" \
    "${fastq_dir}/"*.fastq.gz

# run fastq-screen (threads option inherited by bwa/bowtie)
echo "Running FastQ Screen..."
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

# TODO parallel option {} is incompatible with complex substitutions and variables
# currently requires hard-coding the fastq read suffix
# solution could be to loop through basenames and extend them as necessary, rather than relying on {} syntax
# --out2 '{= s:.*/::; s:\.[^/.]+$::; s:\.[^/.]+$::; s/R1/R2/ =}'.trim.fastq.gz
echo "Running fastp using parallel..."
find "${fastq_dir}" -name *"R1.fastq.gz" |
    parallel -j $((${n_threads} / 2)) --plus \
        fastp \
        --in1 {} \
        --in2 '{= s/_R1/_R2/ =}' \
        --out1 "${output_dir}/fastp/"{/..}.trim.fastq.gz \
        --out2 "${output_dir}/fastp/"'{= s:.*/::; s/_R1.fastq.gz/_R2.trim.fastq.gz/ =}' \
        --json "${output_dir}/fastp/"'{= s:.*/::; s/_R1.fastq.gz/.trim.json/ =}' \
        --html "${output_dir}/fastp/"'{= s:.*/::; s/_R1.fastq.gz/.trim.html/ =}' \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_tail \
        --qualified_quality_phred 20 \
        --length_required 15 \
        --trim_poly_x \
        --thread 2

# re-run qc after trimming
echo "Re-running FastQC after trimming..."
fastqc \
    --threads "${n_threads}" \
    --outdir "${output_dir}/fastqc/" \
    "${output_dir}/fastp/"*".trim.fastq.gz"

# aggregate results with multiQC
multiqc --force "${output_dir}" --config "${multiqc_conf}" --outdir "${output_dir}/multiqc"

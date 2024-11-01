#!/usr/bin/env bash

######################################
# Script to perform variant calling #
######################################

# set bash strict mode
set -euo pipefail

# allow debug mode by running `TRACE=1 ./script.sh` - equivalent to `set -x`
if [[ "${TRACE-0}" == "1" ]]; then set -o xtrace; fi

# get file path of script to allow it to be run from any working directory
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
PROJECT_ROOT=$(realpath "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)/../")
echo "Project root = ${PROJECT_ROOT}"

#####################
# Options and paths #
#####################

# set number of threads for downstream tools
n_threads="${SLURM_CPUS_PER_TASK:-8}"

# define in and outputs
ref_human="${PROJECT_ROOT}/data/ref/human/Homo_sapiens.GRCh38.p14.GENCODE.release45/GRCh38.primary_assembly.genome.fa.gz"
ref_pf="${PROJECT_ROOT}/data/ref/Pfalciparum/PlasmoDB-release-68/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"
ref_pv="${PROJECT_ROOT}/data/ref/Pvivax/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.fasta"
ref_pm="${PROJECT_ROOT}/data/ref/Pmalariae/PlasmoDB-release-68/PlasmoDB-68_PmalariaeUG01_Genome.fasta"
ref_poc="${PROJECT_ROOT}/data/ref/Povale/PlasmoDB-release-68/PlasmoDB-68_PovalecurtisiGH01_Genome.fasta"
ref_pow="${PROJECT_ROOT}/data/ref/Povale/PlasmoDB-release-68/PlasmoDB-68_PovalewallikeriPowCR01_Genome.fasta"

output_dir="${PROJECT_ROOT}/results/"
fastq_dir="${output_dir}/fastp/"
bam_dir="${output_dir}/bwa/"
vcf_dir="${output_dir}/gatk/"
multiqc_conf="${PROJECT_ROOT}/config/multiqc_config.yaml"

# create output directories
mkdir -p "${vcf_dir}"

# check if bam directory exist
if [ ! -d "${bam_dir}" ]; then
    echo "BAM input directory (${bam_dir}) does not exist."
    exit 1
fi

# check if reference fasta files exists
# for ref in ${ref_human} ${ref_pf} ${ref_pv} ${ref_poc} ${ref_pow} ${ref_pm}; do
for ref in ${ref_pf} ${ref_pv} ${ref_poc} ${ref_pow} ${ref_pm}; do
    if ! [ -f "${ref}" ]; then
        echo "Reference fasta file not found (${ref})."
        exit 1
    fi
done

# log run options
printf "
GATK variant calling script | $(basename "${BASH_SOURCE[0]}")
==============================================

Output directory:           ${bam_dir}
FASTQ reads directory:      ${fastq_dir}
Reference human:            ${ref_human}
Reference Pfalciparum:      ${ref_pf}
Reference Pvivax:           ${ref_pv}
Reference Pmalaria:         ${ref_pm}
Reference Povale w:         ${ref_pow}
Refence Povale c            ${ref_poc}
threads:                    ${n_threads}
"

####################
# Start of mapping #
####################

# create reference fai and dict files if they do not yet exist
for ref in ${ref_pf} ${ref_pv} ${ref_pm} ${ref_pow} ${ref_poc}; do
    index_files_found=1
    if ! [ -f "${ref}.fai" ]; then
        index_files_found=0
        echo "Building samtools faidx index for ${ref}..."
        samtools faidx "${ref}"
    fi
    if ! [ -f "${ref%.fasta}.dict" ]; then
        index_files_found=0
        echo "Building GATK sequence dictionary for ${ref}..."
        gatk CreateSequenceDictionary -R "${ref}"
    fi
    if ! [ -f "${ref%.fasta}.bed" ]; then
        index_files_found=0
        echo "Creating region-level bed file for ${ref}..."
        awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${ref}.fai > ${ref%.fasta}.bed
    fi
    if [ "$index_files_found" -eq 1 ]; then
        echo "Found fai, dict and bed files for ${ref}, skipping indexing steps..."
    fi
done

# function to call variants for a single sample
function haplotype_caller() {
    bam=${1}

	sample_name=$(basename ${bam} .sort.markdup.bam)

    species=$(awk -v pat="${sample_name}" -F',' '$1 ~ pat { print $2; exit}' "${PROJECT_ROOT}/data/samplesheet.csv")

    if [[ "${species}" == "pf" ]]; then
        ref="${ref_pf}"
    elif [[ "${species}" == "pv" ]]; then
        ref="${ref_pv}"
    elif [[ "${species}" == "pm" ]]; then
        ref="${ref_pm}"
    elif [[ "${species}" == "pow" ]]; then
        ref="${ref_pow}"
    elif [[ "${species}" == "poc" ]]; then
        ref="${ref_poc}"
    fi

    printf "\nCalling variants for ${bam} (${species})...\n"

    gatk --java-options "-Xmx${mem}g" HaplotypeCaller \
        -R "${ref}" \
        -I "${bam}" \
        -O "${vcf_dir}/${sample_name}.g.vcf.gz" \
        --native-pair-hmm-threads 4 \
        -ERC GVCF
}
export -f haplotype_caller

# parallellize variant calling across samples
jobs=$((${n_threads}/4))
if [ -n "${SLURM_MEM_PER_NODE-}" ]; then
    mem=$((${SLURM_MEM_PER_NODE}/1000/${jobs}))
elif [ -n "${SLURM_MEM_PER_CPU-}" ]; then
    mem=$((${SLURM_MEM_PER_CPU}*4))
else
    mem=4
fi

# make variables required by function available
export PROJECT_ROOT vcf_dir mem ref_pf ref_pk ref_pv ref_pm ref_poc ref_pow

printf "\nParallellizing across ${jobs} jobs and assigning each ${mem}G of memory...\n"

parallel -j "${jobs}" \
    haplotype_caller "{}" \
    ::: "${bam_dir}"/*.sort.markdup.bam

# aggregate results with multiQC
# multiqc --force "${output_dir}" --config "${multiqc_conf}" --outdir "${output_dir}/multiqc"

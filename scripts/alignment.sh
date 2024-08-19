#!/usr/bin/env bash

# set bash strict mode
set -euo pipefail

# allow debug mode by running `TRACE=1 ./script.sh` - equivalent to `set -x`
if [[ "${TRACE-0}" == "1" ]]; then set -o xtrace; fi

# get file path of script
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
# Get file path of script and set project root.
# PROJECT_ROOT=$(realpath "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)/../../")
# echo "Project root = ${PROJECT_ROOT}"

echo "Script dir = ${SCRIPT_DIR}"

#####################
# Options and paths #
#####################

# set number of threads for downstream tools
n_threads="${SLURM_CPUS_PER_TASK:-8}"

# define in and outputs
ref_pf="${SCRIPT_DIR}/../data/ref/Pfalciparum/PlasmoDB-release-68/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"
ref_pv="${SCRIPT_DIR}/../data/ref/Pvivax/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.fasta"
ref_pm="${SCRIPT_DIR}/../data/ref/Pmalariae/PlasmoDB-release-68/PlasmoDB-68_PmalariaeUG01_Genome.fasta"
ref_poc="${SCRIPT_DIR}/../data/ref/Povale/PlasmoDB-release-68/PlasmoDB-68_PovalecurtisiGH01_Genome.fasta"
ref_pow="${SCRIPT_DIR}/../data/ref/Povale/PlasmoDB-release-68/PlasmoDB-68_PovalewallikeriPowCR01_Genome.fasta"
# ref_human="${SCRIPT_DIR}/../data/ref/GRCh38.chr21.fa.gz"
# ref_concat="${SCRIPT_DIR}/../data/ref/concat.fa.gz"

fastq_dir="${SCRIPT_DIR}/../results/fastp/"
# reads_human="${SCRIPT_DIR}/../data/fastq/human"

output_dir="${SCRIPT_DIR}/../results/bwa/"
mkdir -p "${output_dir}"

# check if fastq directory exist
if [ ! -d "${fastq_dir}" ]; then
    echo "FASTQ input directory (${fastq_dir}) does not exist."
fi

# # check if reference fasta file exists
# if ! [ -f "${ref}" ]; then
    # echo "Reference fasta file not found (${ref})."
# fi

# log run options
printf "
BWA MEM script | $(basename "${BASH_SOURCE[0]}")
==============================================

Output directory:           ${output_dir}
FASTQ reads directory:      ${fastq_dir}
Reference Pfalciparum:      ${ref_pf}
Reference Pvivax:           ${ref_pv}
Reference Pmalaria:         ${ref_pm}
Reference Povale w:         ${ref_pow}
Refence Povale c            ${ref_poc}
threads:                    ${n_threads}
"

###################
# start of script #
###################

# create reference index if it does not yet exist
# required for fastq-screen
for ref in ${ref_pf} ${ref_pv} ${ref_pm} ${ref_pow} ${ref_poc}; do
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

# for i in "${ref}."{amb,ann,bwt,pac,sa}; do
    # if ! [ -f "${i}" ]; then
        # bwa index "${ref}"
        # break
    # fi
# done

# map fastq read pairs
for r1 in "${fastq_dir}"/*_R1_001.fastp.fastq.gz; do

    # get filepath containing basename of each read pair
    sample_path="${r1%_R1_001.fastp.fastq.gz}"

    # convert to basename of each read without the filepath prefix
    sample_name="${sample_path##*/}"

    # create output filepath for each read pair
    output_prefix="${output_dir}/${sample_name}"
    # mkdir -p "${output_prefix}"

    printf "\nProcessing sample %s ...\n\n" "${sample_name}"

    species=$(awk -v pat="${sample_name}" -F',' '$1 ~ pat { print $2}' "${SCRIPT_DIR}/samplesheet_bwa.csv")

    if [[ "${species}" == "pf" ]]; then
        ref="${ref_pf}"
    elif [[ "${species}" == "pm" ]]; then
        ref="${ref_pm}"
    elif [[ "${species}" == "pow" ]]; then
        ref="${ref_pow}"
    fi

    echo "Species is $species, ref is $ref"

    echo "Creating bam file: ${output_prefix}.sort.bam"

    # map to reference genome
    bwa mem \
        -t "${n_threads}" \
        -Y -K 100000000 \
        -R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA" \
        "${ref}" \
        "${sample_path}_R1_001.fastp.fastq.gz" \
        "${sample_path}_R2_001.fastp.fastq.gz" |
        # sort and compress to bam
        samtools sort --threads "${n_threads}" \
            -o "${output_prefix}.sort.bam"

    # generate stats
    samtools index --threads "${n_threads}" "${output_prefix}.sort.bam"
    samtools stats --threads "${n_threads}" "${output_prefix}.sort.bam" >"${output_prefix}.sort.bam.stats"
    samtools flagstat --threads "${n_threads}" "${output_prefix}.sort.bam" >"${output_prefix}.sort.bam.flagstat"
    samtools idxstats --threads "${n_threads}" "${output_prefix}.sort.bam" >"${output_prefix}.sort.bam.idxstats"

done

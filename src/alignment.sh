#!/usr/bin/env bash

# set bash strict mode
set -euo pipefail

# allow debug mode by running `TRACE=1 ./script.sh` - equivalent to `set -x`
if [[ "${TRACE-0}" == "1" ]]; then set -o xtrace; fi

# get file path of project root to allow it to be run from any working directory
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
PROJECT_ROOT=$(realpath "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)/../")
echo "Project root = ${PROJECT_ROOT}"

# set number of threads for downstream tools
n_threads=8
if [ -n "${SLURM_MEM_PER_NODE-}" ]; then
    mem_in_gb=$((${SLURM_MEM_PER_NODE}/1000))
elif [ -n "${SLURM_MEM_PER_CPU-}" ]; then
    mem_in_gb=$((${SLURM_MEM_PER_CPU}/1000))
else
    mem_in_gb=8
fi

# define in and outputs
ref="${PROJECT_DIR}/data/ref/Pf3D7_01_v3.fa"
output_dir="${PROJECT_ROOT}/results/"
trimmed_fastq_dir="${output_dir}/fastp/"
bam_dir="${output_dir}/bwa/"
mkdir -p "${bam_dir}"

# check if fastq directory exist
if [ ! -d "${trimmed_fastq_dir}" ]; then
    echo "Trimmed FASTQ directory (${trimmed_fastq_dir}) does not exist."
fi

# check if reference fasta file exists
if ! [ -f "${ref}" ]; then
    echo "Reference fasta file not found (${ref})."
fi

# log run options
printf "
BWA MEM script | $(basename "${BASH_SOURCE[0]}")
==============================================

Output directory:           ${bam_dir}
Trimmed FASTQ reads:        ${trimmed_fastq_dir}
Reference pfalciparum:      ${ref}
threads:                    ${n_threads}
"

####################
# Start of mapping #
####################

# create reference index if it does not yet exist
for i in "${ref}."{amb,ann,bwt,pac,sa}; do
    if ! [ -f "${i}" ]; then
        echo "Building BWA index for ${ref}..."
        bwa index "${ref}"
        break
    fi
done

# map fastq read pairs using bwa
for r1 in "${trimmed_fastq_dir}"/*_R1_001.trim.fastq.gz; do

    # get filepath containing basename of each read pair
    sample_path="${r1%_R1_001.trim.fastq.gz}"

    # convert to basename of each read without the filepath prefix
    sample_name="${sample_path##*/}"

    # create output filepath for each read pair
    output_prefix="${bam_dir}/${sample_name}"

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

# Picard mark duplicates
printf "\nMarking duplicates...\n"

for bam in "${bam_dir}/"*.sort.bam; do
    sample_name=$(basename "${bam}" ".sort.bam")
    gatk --java-options "-Xmx${mem_in_gb}g" \
        MarkDuplicates \
        --INPUT "${bam}" \
        --OUTPUT "${bam_dir}/${sample_name}.sort.markdup.bam" \
        --METRICS_FILE "${bam_dir}/${sample_name}.sort.markdup.metrics" \
        --REMOVE_DUPLICATES false
done

# Create summary statistics
for bam in "${bam_dir}/"*.sort.markdup.bam; do
    printf "\nCreating index and samtool stats for ${bam}...\n"
    samtools index --threads "${n_threads}" "${bam}"
    samtools stats --threads "${n_threads}" "${bam}" >"${bam}.stats"
    samtools flagstat --threads "${n_threads}" "${bam}" >"${bam}.flagstat"
    samtools idxstats --threads "${n_threads}" "${bam}" >"${bam}.idxstats"
done

# clean up
# rm "${bam_dir}/"*.sort.bam

# aggregate results with multiQC
multiqc --force "${output_dir}" --config "${multiqc_conf}" --outdir "${output_dir}/multiqc"

#!/usr/bin/env bash

#################################################################
# Script to perform quality control and trimming of fastq reads #
#################################################################

# set bash strict mode
set -euo pipefail

# allow debug mode by running `TRACE=1 ./script.sh` - equivalent to `set -x`
if [[ "${TRACE-0}" == "1" ]]; then set -o xtrace; fi

# get file path of project root to allow it to be run from any working directory
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
PROJECT_ROOT=$(realpath "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)/../")
echo "Project root = ${PROJECT_ROOT}"

#####################
# Options and paths #
#####################

# set number of threads and memory for downstream tools
n_threads="${SLURM_CPUS_PER_TASK:-8}"

# define in and outputs
output_dir="${PROJECT_ROOT}/results/"
trimmed_fastq_dir="${output_dir}/fastp/"
bam_dir="${output_dir}/bwa/"

# create output directories
mkdir -p "${bam_dir}"

# read config files
multiqc_conf="${PROJECT_ROOT}/config/multiqc_config.yaml"

# TODO: define fastq read suffix
# read_1_suffix="_R1_001.fastq.gz"
# read_2_suffix="_R2_001.fastq.gz"
# read_file_extension=".fastq.gz"

# reference files
# TODO: either make symlinks from the reference genome location to the
# data/ref/ directory in the project: e.g.,
# `ln -s /data/antwerpen/grp/ap_itg_mu/public_data/reference_genomes/Pvivax/PlasmoDB-release-68/* ./data/ref/`
# or use the full path to the reference genomes
# ref="${PROJECT_ROOT}/data/ref/GCF_000001405.40_GRCh38.p14_genomic-PlasmoDB-68_PvivaxPAM_Genome.fa.bgz"
# ref="${PROJECT_ROOT}/data/ref/GENCODE.release45-GRCh38.primary_assembly.genome-PlasmoDB-68_PvivaxPAM_Genome.fa.bgz"
# ref="/data/antwerpen/grp/ap_itg_mu/public_data/reference_genomes/combined/human_pv/GCF_000001405.40_GRCh38.p14_genomic-PlasmoDB-68_PvivaxPAM_Genome.fa.bgz"
ref="/data/antwerpen/grp/ap_itg_mu/public_data/reference_genomes/combined/human_pv/GENCODE.release45-GRCh38.primary_assembly.genome-PlasmoDB-68_PvivaxPAM_Genome.fa.bgz"
# ref_pv="${PROJECT_ROOT}/data/ref/Pvivax/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.fasta"
ref_pv="/data/antwerpen/grp/ap_itg_mu/public_data/reference_genomes/Pvivax/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.fasta"
# TODO: what to do if file extension is not fasta but fa or fa.gz?
bed_pv="${ref_pv%.fasta}.bed"
# bed_pv="${ref_pv}.bed"

# check if trimmed fastq files exist
if [ ! -d "${trimmed_fastq_dir}" ]; then
    printf "\nTrimmed FASTQ directory (${trimmed_fastq_dir}) does not exist.\n"
    exit 1
fi

# check if reference fasta files exists
for i in ${ref} ${ref_pv}; do
    if ! [ -f "${i}" ]; then
        printf "\nReference fasta file not found (${i}).\n"
        exit 1
    fi
done

# log run options
printf "
BWA MEM script | $(basename "${BASH_SOURCE[0]}")
==============================================

Output directory:           ${bam_dir}
Trimmed FASTQ reads:        ${trimmed_fastq_dir}
Reference:                  ${ref}
threads:                    ${n_threads}
"

####################
# Start of mapping #
####################

# create reference index if it does not yet exist
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

# Create bed file with chromsomes/regions based on Pv fasta file
if ! [ -f "${bed_pv}" ]; then
    samtools faidx "${ref_pv}"
    awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' "${ref_pv}.fai" > "${bed_pv}"
else
    printf "\nFound .bed file for ${ref_pv}, skipping .bed file generation step...\n"
fi

# map fastq read pairs using bwa - competitive against concatenated human-plasmodium reference
for r1 in "${trimmed_fastq_dir}"/*_R1.trim.fastq.gz; do

    # get filepath containing basename of each read pair
    sample_path="${r1%_R1.trim.fastq.gz}"

    # convert to basename of each read without the filepath prefix
    sample_name="${sample_path##*/}"

    # TODO: modify this to read from samplesheet.csv instead
    # retrieve lane and sample group
    sample_lane=$(echo ${sample_name} | grep -Po 'L\d{3}')
    sample=$(echo ${sample_name} | cut -d '_' -f2)
    sample_flowcell=$(echo ${sample_name} | cut -d '_' -f1)
    sample_library=${sample}

    printf "\nMapping raw reads to combined human-plasmodium reference for sample %s, lane %s of sample group %s ...\n\n" "${sample_name}" "${sample_lane}" "${sample_library}"
    printf "\nCreating combined bam file: %s.sort.human.bam\n" "${bam_dir}/${sample_name}"

    bwa mem \
        -t "${n_threads}" \
        -Y -K 100000000 \
        -R "@RG\tID:${sample_name}\tSM:${sample}\tPL:ILLUMINA\\tPU:${sample_flowcell}.${sample_lane}\\tLB:${sample_library}" \
        "${ref}" \
        "${sample_path}_R1.trim.fastq.gz" \
        "${sample_path}_R2.trim.fastq.gz" |
    # sort and compress to bam
        samtools sort --threads "${n_threads}" \
            -o "${bam_dir}/${sample_name}.sort.combined.bam"

    # extract plasmodium reads
    samtools view -b -h -L "${bed_pv}" "${bam_dir}/${sample_name}.sort.combined.bam" > "${bam_dir}/${sample_name}.sort.plasmodium.bam"
done

# picard mark duplicates
jobs=$((${n_threads}/2))
if [ -n "${SLURM_MEM_PER_NODE-}" ]; then
    mem=$((${SLURM_MEM_PER_NODE}/1000/${jobs}))
elif [ -n "${SLURM_MEM_PER_CPU-}" ]; then
    mem=$((${SLURM_MEM_PER_CPU}*2/1000))
else
    mem=4
fi

printf "\nMarking duplicates in parallel using %s GB of memory per job (n_jobs = %s)\n" "${mem}" "${jobs}"

# mark duplicates simultaneously on same sample run across different lanes and merge bam files
function add_input_prefix() {
    declare -a arr=()
    for i in "${2}/"*"${1}"*.sort.plasmodium.bam; do
        arr+=( "--INPUT ${i}" )
    done;
    echo ${arr[@]}
}
export -f add_input_prefix
export bam_dir

# ! Note that single quotes are required to avoid the command substitution around the
# ! add_input_prefix function call from being executed before it is passed to parallel
for i in "${bam_dir}/"*.sort.plasmodium.bam; do echo $(basename $i | cut -d '_' -f2); done | \
    sort -u | \
    # while read -r line; do add_input_prefix $line "${bam_dir}/"; |
    parallel -j "${jobs}" \
        gatk --java-options -Xmx${mem}G \
            MarkDuplicates \
            '$(add_input_prefix {} ${bam_dir})' \
            --OUTPUT "${bam_dir}/{}.sort.plasmodium.markdup.bam" \
            --METRICS_FILE "${bam_dir}/{}.sort.plasmodium.markdup.metrics" \
            --REMOVE_DUPLICATES false

for bam in "${bam_dir}/"*.sort.plasmodium.markdup.bam; do
    printf "\nCreating index and samtool stats for ${bam}...\n"
    samtools index --threads "${n_threads}" "${bam}"
    samtools stats --threads "${n_threads}" "${bam}" >"${bam}.stats"
    samtools flagstat --threads "${n_threads}" "${bam}" >"${bam}.flagstat"
    samtools idxstats --threads "${n_threads}" "${bam}" >"${bam}.idxstats"
done

# clean up
# rm "${bam_dir}/"*.sort.combined.bam

# aggregate results with multiQC
multiqc --force "${output_dir}" --config "${multiqc_conf}" --outdir "${output_dir}/multiqc"

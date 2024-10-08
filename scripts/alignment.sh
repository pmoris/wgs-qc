#!/usr/bin/env bash

#################################################################
# Script to perform quality control and trimming of fastq reads #
#################################################################

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
multiqc_conf="${PROJECT_ROOT}/config/multiqc_config.yaml"
# reads_human="${PROJECT_ROOT}/data/fastq/human"

# create output directories
mkdir -p "${bam_dir}"

# check if fastq directory exist
if [ ! -d "${fastq_dir}" ]; then
    echo "FASTQ input directory (${fastq_dir}) does not exist."
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
BWA MEM script | $(basename "${BASH_SOURCE[0]}")
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

# create reference index if it does not yet exist
# required for fastq-screen
# for ref in ${ref_human} ${ref_pf} ${ref_pv} ${ref_pm} ${ref_pow} ${ref_poc}; do
for ref in ${ref_human} ${ref_pf} ${ref_pv} ${ref_pm} ${ref_pow} ${ref_poc}; do
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

# map fastq read pairs using bwa
for r1 in "${fastq_dir}"/*_R1_001.trim.fastq.gz; do

    # get filepath containing basename of each read pair
    sample_path="${r1%_R1_001.trim.fastq.gz}"

    # convert to basename of each read without the filepath prefix
    sample_name="${sample_path##*/}"

    # TODO: modify this to read from samplesheet.csv instead
    # retrieve lane and sample group
    sample_lane="${sample_name##*_}"
    sample_group="${sample_name%%_L*}"
    sample_flowcell="${sample_group##*_}"
    sample="${sample_name%%_*}"

    # create output filepath for each read pair
    output_prefix="${bam_dir}/${sample_name}"

    species=$(awk -v pat="${sample_name}" -F',' '$1 ~ pat { print $2}' "${PROJECT_ROOT}/data/samplesheet.csv")

    if [[ "${species}" == "pf" ]]; then
        ref="${ref_pf}"
    elif [[ "${species}" == "pv" ]]; then
        ref="${ref_pv}"
    elif [[ "${species}" == "pm" ]]; then
        ref="${ref_pm}"
    elif [[ "${species}" == "pow" ]]; then
        ref="${ref_pow}"
    fi

    # map to human reference genome first to remove host reads
    printf "\nMapping raw reads to human reference for sample %s, lane %s of sample group %s ...\n\n" "${sample_name}" "${sample_lane}" "${sample_group}"
    echo "Creating human bam file: ${output_prefix}.sort.human.bam"
    bwa mem \
        -t "${n_threads}" \
        -Y -K 100000000 \
        -R "@RG\tID:${sample_name}\tSM:${sample}\tPL:ILLUMINA\\tPU:${sample_flowcell}.${sample_lane}\\tLB:${sample_group}" \
        "${ref_human}" \
        "${sample_path}_R1_001.trim.fastq.gz" \
        "${sample_path}_R2_001.trim.fastq.gz" |
    # sort and compress to bam
        samtools sort --threads "${n_threads}" \
            -o "${output_prefix}.sort.human.bam"

    # extract all unmapped pairs (both reads unmapped)
    # approach adapted from https://lh3.github.io/2021/07/06/remapping-an-aligned-bam
    # TODO: alternatively use bedtools' bamtofastq approach and save intermediate steps
    printf "\nMapping filtered reads to %s genome for sample %s, lane %s of sample group %s ...\n\n" "${species}" "${sample_name}" "${sample_lane}" "${sample_group}"
    echo "Creating bam file: ${output_prefix}.sort.bam"
    samtools view -b -f 12 "${output_prefix}.sort.human.bam" |
    # convert back to fastq
        samtools collate -Oun128 - |
        samtools fastq -OT RG,BC - |
    # map to plasmodium genome
    # -CH adds back original read group info
    # -p gathers paired reads from stream - https://github.com/samtools/samtools/issues/1306
        bwa mem \
            -t "${n_threads}" \
            -Y -K 100000000 \
            -CH <(samtools view -H "${output_prefix}.sort.human.bam" | grep ^@RG) \
            -p \
            "${ref}" \
            - |
    # sort and compress to bam
        samtools sort --threads "${n_threads}" \
            -o "${output_prefix}.sort.bam"

    # simple mapping method without host removal
    # # map to reference genome
    # bwa mem \
    #     -t "${n_threads}" \
    #     -Y -K 100000000 \
    #     -R "@RG\tID:${sample_name}\tSM:${sample}\tPL:ILLUMINA\\tPU:${sample_flowcell}.${sample_lane}\\tLB:${sample_group}" \
    #     "${ref}" \
    #     "${sample_path}_R1_001.trim.fastq.gz" \
    #     "${sample_path}_R2_001.trim.fastq.gz" |
    # # sort and compress to bam
    #     samtools sort --threads "${n_threads}" \
    #         -o "${output_prefix}.sort.bam"

    # # generate stats
    # samtools index --threads "${n_threads}" "${output_prefix}.sort.bam"
    # samtools stats --threads "${n_threads}" "${output_prefix}.sort.bam" >"${output_prefix}.sort.bam.stats"
    # samtools flagstat --threads "${n_threads}" "${output_prefix}.sort.bam" >"${output_prefix}.sort.bam.flagstat"
    # samtools idxstats --threads "${n_threads}" "${output_prefix}.sort.bam" >"${output_prefix}.sort.bam.idxstats"

done

# picard mark duplicates
printf "\nMarking duplicates...\n"

jobs=$((${n_threads}/2))

# mark duplicates for each file separately -> requires manual merging afterwards
# parallel -j "${jobs}" \
#     gatk --java-options -Xmx$((8))G \
#         MarkDuplicates \
#         --INPUT "{}" \
#         --OUTPUT "${bam_dir}/{/.}.markdup.bam" \
#         --METRICS_FILE {.}.markdup.metrics \
#         --REMOVE_DUPLICATES false \
#     ::: "${bam_dir}"/*.sort.bam

# mark duplicates simultaneously on same sample run across different lanes with automatic concatenation

# ! Note that the internal loop of the prefix function needs to
# ! use a glob pattern that includes _L, to avoid libraries sharing
# ! a prefix ID from being grouped together. E.g., sample_S11_L001 and sample S_1_L001.
# --INPUT fastq/ANT5670_S11_L001_R1_001.fastq.gz --INPUT fastq/ANT5670_S11_L001_R2_001.fastq.gz --INPUT fastq/ANT5670_S1_L001_R1_001.fastq.gz --INPUT fastq/ANT5670_S1_L001_R2_001.fastq.gz

# testing:
# $ for i in "fastq/"*.fastq.gz; do echo "${i%%_L*}"; done | sort -u | while read -r line; do declare -a arr=(); for i in $line*; do arr+=( "--INPUT ${i}" ); done; echo ${arr[@]}; done
# --INPUT fastq/ANT5670_S11_L001_R1_001.fastq.gz --INPUT fastq/ANT5670_S11_L001_R2_001.fastq.gz --INPUT fastq/ANT5670_S1_L001_R1_001.fastq.gz --INPUT fastq/ANT5670_S1_L001_R2_001.fastq.gz
# --INPUT fastq/ANT5670_S11_L001_R1_001.fastq.gz --INPUT fastq/ANT5670_S11_L001_R2_001.fastq.gz
# --INPUT fastq/ANT5670_S261_L001_R1_001.fastq.gz --INPUT fastq/ANT5670_S261_L001_R2_001.fastq.gz --INPUT fastq/ANT5670_S261_L002_R1_001.fastq.gz --INPUT fastq/ANT5670_S261_L002_R2_001.fastq.gz
# --INPUT fastq/ANT6000_S1_L001_R1_001.fastq.gz --INPUT fastq/ANT6000_S1_L001_R2_001.fastq.gz

# $ for i in "fastq/"*.fastq.gz; do echo "${i%%_L*}"; done | sort -u | while read -r line; do declare -a arr=(); for i in ${line}_L*; do arr+=( "--INPUT ${i}" ); done; echo ${arr[@]}; done
# --INPUT fastq/ANT5670_S1_L001_R1_001.fastq.gz --INPUT fastq/ANT5670_S1_L001_R2_001.fastq.gz
# --INPUT fastq/ANT5670_S11_L001_R1_001.fastq.gz --INPUT fastq/ANT5670_S11_L001_R2_001.fastq.gz
# --INPUT fastq/ANT5670_S261_L001_R1_001.fastq.gz --INPUT fastq/ANT5670_S261_L001_R2_001.fastq.gz --INPUT fastq/ANT5670_S261_L002_R1_001.fastq.gz --INPUT fastq/ANT5670_S261_L002_R2_001.fastq.gz
# --INPUT fastq/ANT6000_S1_L001_R1_001.fastq.gz --INPUT fastq/ANT6000_S1_L001_R2_001.fastq.gz

# alternatively, change initial loop to report including _ before L (and also L itself?) and then remove it again when creating outputs?

function add_input_prefix() {
    # echo ${1};
    declare -a arr=()
    for i in "${1}_L"*".sort.bam"; do
        # echo "looping over ${i}"
        arr+=( "--INPUT ${i}" )
    done;
    echo ${arr[@]}
}
export -f add_input_prefix

# ! Note that single quotes are required to avoid the command substitution around the
# ! add_input_prefix function call from being executed before it is passed to parallel
for i in "${bam_dir}/"*.sort.bam; do echo "${i%%_L*}"; done | sort -u | \
    parallel -j "${jobs}" \
        gatk --java-options -Xmx$((8))G \
            MarkDuplicates \
            '$(add_input_prefix {})' \
            --OUTPUT "{}.sort.markdup.bam" \
            --METRICS_FILE {}.markdup.metrics \
            --REMOVE_DUPLICATES false

# for bam in results-testset/bwa/*.sort.bam; do echo "${bam%%_L*}"; done | sort -u | while read -r line ; do array=(${line}*.sort.bam); echo ${array[@]}; done

# └─▶ for bam in results-testset/bwa/*.sort.bam; do echo "${bam%%_L*}"; done | sort -u | while read -r line ; do echo $(cat "sdfgsd" ${line}*.sort.bam); done

# └─▶ for i in results-testset/bwa/*.sort.bam; do echo "${i%%_L*}"; done | sort -u | parallel echo "{}*.sort.bam" "{}"

# # TODO: incorporate into parallel?
for bam in "${bam_dir}/"*.sort.markdup.bam; do
    printf "\nCreating index and samtool stats for ${bam}...\n"
    samtools index --threads "${n_threads}" "${bam}"
    samtools stats --threads "${n_threads}" "${bam}" >"${bam}.stats"
    samtools flagstat --threads "${n_threads}" "${bam}" >"${bam}.flagstat"
    samtools idxstats --threads "${n_threads}" "${bam}" >"${bam}.idxstats"
done

for bam in "${bam_dir}/"*.sort.human.bam; do
    printf "\nCreating index and samtool stats for ${bam}...\n"
    samtools index --threads "${n_threads}" "${bam}"
    samtools stats --threads "${n_threads}" "${bam}" >"${bam}.stats"
    samtools flagstat --threads "${n_threads}" "${bam}" >"${bam}.flagstat"
    samtools idxstats --threads "${n_threads}" "${bam}" >"${bam}.idxstats"
done

# clean up
rm "${bam_dir}/"*.sort.bam "${bam_dir}/"*.sort.human.bam

# aggregate results with multiQC
multiqc --force "${output_dir}" --config "${multiqc_conf}" --outdir "${output_dir}/multiqc"

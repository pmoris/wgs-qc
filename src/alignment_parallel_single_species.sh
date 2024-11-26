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

# read species from command line
species=$1

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

# check if trimmed fastq files exist
if [ ! -d "${trimmed_fastq_dir}" ]; then
    echo "Trimmed FASTQ directory (${trimmed_fastq_dir}) does not exist."
    exit 1
fi

# check if reference fasta files exists
for ref in ${ref_human} ${ref_plasmodium}; do
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
Trimmed FASTQ reads:        ${trimmed_fastq_dir}
Reference human:            ${ref_human}
Reference Plasmodium:           ${ref_plasmodium}
threads:                    ${n_threads}
"

####################
# Start of mapping #
####################

# create reference index if it does not yet exist
for ref in ${ref_human} ${ref_plasmodium}; do
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

    # map to human reference genome first to remove host reads
    printf "\nMapping raw reads to human reference for sample %s, lane %s of sample group %s ...\n\n" "${sample_name}" "${sample_lane}" "${sample_library}"
    printf "Creating human bam file: %s.sort.human.bam" "${bam_dir}/${sample_name}"
    bwa mem \
        -t "${n_threads}" \
        -Y -K 100000000 \
        -R "@RG\tID:${sample_name}\tSM:${sample}\tPL:ILLUMINA\\tPU:${sample_flowcell}.${sample_lane}\\tLB:${sample_library}" \
        "${ref_human}" \
        "${sample_path}_R1.trim.fastq.gz" \
        "${sample_path}_R2.trim.fastq.gz" |
    # sort and compress to bam
        samtools sort --threads "${n_threads}" \
            -o "${bam_dir}/${sample_name}.sort.human.bam"

    # extract all unmapped pairs (both reads unmapped)
    # approach adapted from https://lh3.github.io/2021/07/06/remapping-an-aligned-bam
    # TODO: alternatively use bedtools' bamtofastq approach and save intermediate steps
    printf "\nMapping filtered reads to Plasmodium genome for sample %s, lane %s of sample group %s ...\n\n" "${sample_name}" "${sample_lane}" "${sample_library}"
    printf "Creating bam file: %s.sort.bam" "${bam_dir}/${sample_name}"
    samtools view -b -f 12 "${bam_dir}/${sample_name}.sort.human.bam" |
    # convert back to fastq
        samtools collate -Oun128 - |
        samtools fastq -OT RG,BC - |
    # map to plasmodium genome
    # -CH adds back original read group info
    # -p gathers paired reads from stream - https://github.com/samtools/samtools/issues/1306
        bwa mem \
            -t "${n_threads}" \
            -Y -K 100000000 \
            -CH <(samtools view -H "${bam_dir}/${sample_name}.sort.human.bam" | grep ^@RG) \
            -p \
            "${ref_plasmodium}" \
            - |
    # sort and compress to bam
        samtools sort --threads "${n_threads}" \
            -o "${bam_dir}/${sample_name}.sort.bam"
done

# picard mark duplicates
printf "\nMarking duplicates...\n"

jobs=$((${n_threads}/2))
if [ -n "${SLURM_MEM_PER_NODE-}" ]; then
    mem=$((${SLURM_MEM_PER_NODE}/1000/${jobs}))
elif [ -n "${SLURM_MEM_PER_CPU-}" ]; then
    mem=$((${SLURM_MEM_PER_CPU}*2/1000))
else
    mem=4
fi

printf "\nUsing %s GB of memory per job (n_jobs = %s)\n" "${mem}" "${jobs}"

# mark duplicates simultaneously on same sample run across different lanes and merge bam files
function add_input_prefix() {
    declare -a arr=()
    for i in "${2}/"*"${1}"*.sort.bam; do
        arr+=( "--INPUT ${i}" )
    done;
    echo ${arr[@]}
}
export -f add_input_prefix
export bam_dir

# ! Note that single quotes are required to avoid the command substitution around the
# ! add_input_prefix function call from being executed before it is passed to parallel
echo ${bam_dir}
for i in "${bam_dir}/"*.sort.bam; do echo $(basename $i | cut -d '_' -f2); done | \
    sort -u | \
    # while read -r line; do add_input_prefix $line "${bam_dir}/"; |
    parallel -j "${jobs}" \
        gatk --java-options -Xmx${mem}G \
            MarkDuplicates \
            '$(add_input_prefix {} ${bam_dir})' \
            --OUTPUT "${bam_dir}/{}.sort.markdup.bam" \
            --METRICS_FILE "${bam_dir}/{}.markdup.metrics" \
            --REMOVE_DUPLICATES false

# generate stats
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
# rm "${bam_dir}/"*.sort.bam "${bam_dir}/"*.sort.human.bam

# aggregate results with multiQC
multiqc --force "${output_dir}" --config "${multiqc_conf}" --outdir "${output_dir}/multiqc"

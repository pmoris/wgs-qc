#!/usr/bin/env bash

#################################################################
# Script to perform quality control and trimming of fastq reads #
#################################################################

# set bash strict mode
set -euo pipefail
# set -eo pipefail

# allow debug mode by running `TRACE=1 ./script.sh` - equivalent to `set -x`
if [[ "${TRACE-0}" == "1" ]]; then set -o xtrace; fi

# get file path of project root to allow it to be run from any working directory
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
PROJECT_ROOT=$(realpath "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)/../")
echo "Project root = ${PROJECT_ROOT}"

#####################
# Options and paths #
#####################

die() {
    printf '%s\n' "$1" >&2
    exit 1
}

show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-s SAMPLESHEET.CSV ] [-o OUTPUT DIRECTORY ]
                [-r1 READ 1 SUFFIX ] [-r2 READ 2 SUFFIX ] [-e READ FILE EXTENSION ]
                [-n <name_first/flowcell_first ]
    -h                                      display this help and exit
    -s | --samplesheet SAMPLESHEET.CSV     File path to samplesheet with sample-species info
    -o | --output_dir OUTPUT DIRECTORY      File path to output directory; should already
                                            contain trimmed reads directory named fastp
                                            (default = PROJECT_ROOT/results/)
    -r1 | --read_1_suffix R1_001            Suffix for read pair 1 (excluding file extension)
    -r2 | --read_2_suffix R2_001            Suffix for read pair 2 (excluding file extension)
    -e | --read_file_extension .fastq.gz    Read file extension
    -n | --fastq_identifier <string>        Specifies structure of fastq file name. Either
                                            "name_first" or "flowcell_first".
EOF
}

while :; do
    case ${1:-} in
        -h|-\?|--help)
            show_help    # Display a usage synopsis.
            exit
            ;;

        -s|--samplesheet)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                samplesheet=$2
                shift
            else
                die 'ERROR: "--samplesheet" requires a non-empty option argument.'
            fi
            ;;
        --samplesheet=?*)
            samplesheet=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --samplesheet=)         # Handle the case of an empty --samplesheet=
            die 'ERROR: "--samplesheet" requires a non-empty option argument.'
            ;;


        -o|--output_dir)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                output_dir=$2
                shift
            else
                die 'ERROR: "--output_dir" requires a non-empty option argument.'
            fi
            ;;
        --output_dir=?*)
            output_dir=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --output_dir=)         # Handle the case of an empty --output_dir=
            die 'ERROR: "--output_dir" requires a non-empty option argument.'
            ;;

        -r1|--read_1_suffix)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                read_1_suffix=$2
                shift
            else
                die 'ERROR: "--read_1_suffix" requires a non-empty option argument.'
            fi
            ;;
        --read_1_suffix=?*)
            read_1_suffix=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --read_1_suffix=)         # Handle the case of an empty --output_dir=
            die 'ERROR: "--read_1_suffix" requires a non-empty option argument.'
            ;;

        -r2|--read_2_suffix)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                read_2_suffix=$2
                shift
            else
                die 'ERROR: "--read_2_suffix" requires a non-empty option argument.'
            fi
            ;;
        --read_2_suffix=?*)
            read_2_suffix=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --read_2_suffix=)         # Handle the case of an empty --output_dir=
            die 'ERROR: "--read_2_suffix" requires a non-empty option argument.'
            ;;

        -e|--read_file_extension)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                read_file_extension=$2
                shift
            else
                die 'ERROR: "--read_file_extension" requires a non-empty option argument.'
            fi
            ;;
        --read_file_extension=?*)
            read_file_extension=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --read_file_extension=)         # Handle the case of an empty --output_dir=
            die 'ERROR: "--read_file_extension" requires a non-empty option argument.'
            ;;

        -n|--fastq_identifier)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                fastq_identifier=$2
                shift
            else
                die 'ERROR: "--fastq_identifier" requires a non-empty option argument.'
            fi
            ;;
        --fastq_identifier=?*)
            fastq_identifier=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --fastq_identifier=)         # Handle the case of an empty --output_dir=
            die 'ERROR: "--fastq_identifier" requires a non-empty option argument.'
            ;;

        --)              # End of all options.
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            ;;
        *)               # Default case: No more options, so break out of the loop.
            break
    esac

    shift
done

# set number of threads and memory for downstream tools
n_threads="${SLURM_CPUS_PER_TASK:-8}"

# define default in and outputs
# by defining everything here, there is no need to `cd` to directories first
# ! trailing slash is needed for `find`
samplesheet="$(realpath "${samplesheet:-"${PROJECT_ROOT}/data/samplesheet.csv"}")"
output_dir="$(realpath "${output_dir:-"${PROJECT_ROOT}/results/"}")"
trimmed_fastq_dir="${output_dir}/fastp/"
bam_dir="${output_dir}/bwa/"

# create output directories
mkdir -p "${bam_dir}"

# define default pair suffix and file extensions
read_1_suffix=${read_1_suffix:-"_R1_001"}
read_2_suffix=${read_2_suffix:-"_R2_001"}
read_file_extension=${read_file_extension:-".fastq.gz"}     # extension of trimmed reads should already be set to .fastq.gz in previous script

# set fastq identifier structure
fastq_identifier=${fastq_identifier:-}
if ! [[ "${fastq_identifier}" == "name_first" || "${fastq_identifier}" == "flowcell_first" ]]; then
    printf "\nFastq identifier structure was not set correctly, please specify "name_first" or "flowcell_first" \n"
    exit 1
fi

# config files
multiqc_conf="${PROJECT_ROOT}/config/multiqc_config.yaml"

# reference files
ref_human="${PROJECT_ROOT}/data/ref/human/GRCh38.p14/GENCODE/47/GRCh38.primary_assembly.genome.fa.gz"
ref_pf="${PROJECT_ROOT}/data/ref/Pfalciparum/3D7/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"
ref_pv="${PROJECT_ROOT}/data/ref/Pvivax/PvPAM/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.fasta"
ref_pm="${PROJECT_ROOT}/data/ref/Pmalariae/UG01/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PmalariaeUG01_Genome.fasta"
ref_poc="${PROJECT_ROOT}/data/ref/Povale/curtisi/PocGH01/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PovalecurtisiGH01_Genome.fasta"
ref_pow="${PROJECT_ROOT}/data/ref/Povale/wallikeri/PowCR01/PlasmoDB-release-68/PlasmoDB-68_PovalewallikeriPowCR01_Genome.fasta"
#ref_pk="${PROJECT_ROOT}/data/ref/Pknowlesi/H/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PknowlesiH_Genome.fasta"
ref_phix="${PROJECT_ROOT}/data/ref/PhiX/PhiX-NC_001422.1.fasta"

# check if samplesheet exist
if [ ! -f "${samplesheet}" ]; then
    printf "\nSamplesheet (${samplesheet}) does not exist.\n"
    exit 1
fi

# check if trimmed fastq files exist
if [ ! -d "${trimmed_fastq_dir}" ]; then
    printf "\nTrimmed FASTQ directory (${trimmed_fastq_dir}) does not exist.\n"
    exit 1
fi

# check if reference fasta files exists
for ref in ${ref_human} ${ref_pf} ${ref_pv} ${ref_poc} ${ref_pow} ${ref_pm}; do
    if ! [ -f "${ref}" ]; then
        printf "\nReference fasta file not found (${ref}).\n"
        exit 1
    fi
done

# log run options
printf "
BWA MEM script | $(basename "${BASH_SOURCE[0]}")
==============================================

Output directory:           ${bam_dir}
Samplesheet:                ${samplesheet}
Trimmed FASTQ reads:        ${trimmed_fastq_dir}
Read 1 suffix:              ${read_1_suffix}
Read 2 suffix:              ${read_2_suffix}
Read file extension:        ${read_file_extension}
Fastq identifier:           ${fastq_identifier}
Reference human:            ${ref_human}
Reference Pfalciparum:      ${ref_pf}
Reference Pvivax:           ${ref_pv}
Reference Pmalaria:         ${ref_pm}
Reference Povale wallikeri: ${ref_pow}
Reference Povale curtisi:   ${ref_poc}
threads:                    ${n_threads}
"

####################
# Start of mapping #
####################

# create reference index if it does not yet exist
for ref in ${ref_human} ${ref_pf} ${ref_pv} ${ref_pm} ${ref_pow} ${ref_poc}; do
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


# map fastq read pairs using bwa
for r1 in "${trimmed_fastq_dir}"/*${read_1_suffix}.trim.fastq.gz; do

    # get filepath containing basename of each read pair
    read_file_path="${r1%${read_1_suffix}.trim.fastq.gz}"
    # read_file_name=$(basename "${read_1}" "_R1_001.trim.fastq.gz")

    # convert to basename of each read without the filepath prefix
    read_file_name="${read_file_path##*/}"

    # TODO: modify this to read from samplesheet.csv instead
    # retrieve lane and sample group

    if [[ "${fastq_identifier}" == "name_first" ]]; then

        # ANT5797_S262_L001_R1_001.fastq
        sample_name="$(echo ${read_file_name} | cut -d '_' -f1)"
        # sample="${read_file_name%%_*}"
        sample_lane="$(echo ${read_file_name} | grep -Po 'L\d{3}')"
        # sample_lane="${read_file_name##*_}"
        sample_flowcell="$(zcat "${r1}" | head -n 1 | cut -d ':' -f3)" || true
        sample_library="${sample_name}"

        # it is unclear whether or not the S### identifier refers to unique libraries or not
        # sample_group="${read_file_name%%_L*}"
        # sample_library="${sample_group##*_}"


    elif [[ "${fastq_identifier}" == "flowcell_first" ]]; then

        # 22NY35LT3_106264-002-098_CCTCCTTT-CTTTCGCG_L007_R2.fastq.gz
        sample_name="$(echo ${read_file_name} | cut -d '_' -f2)"
        sample_lane="$(echo ${read_file_name} | grep -Po 'L\d{3}')"
        sample_flowcell="$(echo ${read_file_name} | cut -d '_' -f1)"
        sample_library="${sample_name}"

    fi

    RG_ID="${sample_name}.${sample_flowcell}.${sample_lane}" #.barcode?
    RG_SM="${sample_name}"
    RG_PU="${sample_flowcell}.${sample_lane}"
    RG_LB="${sample_name}.${sample_library}"

    # unset species to make sure there are no left overs from previous iterations
    species=
    species=$(awk -v pat="${sample_name}" -F',' '$1 ~ pat { print $2; exit}' "${samplesheet}")
    if [ -z "${species:-}" ]; then
        echo "Could not find sample during species lookup in samplesheet ${samplesheet}. Exiting..."
        exit 1
    fi

    # unset ref to make sure there are no left overs from previous ref loop
    ref=
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
    if [ -z "${ref:-}" ]; then
        echo "Could not find correct reference based on species lookup in samplesheet ${samplesheet}. Exiting..."
        exit 1
    fi

    # map to human reference genome first to remove host reads
    printf "\nMapping raw reads to human reference for read file %s, flowcell %s, lane %s of sample %s, assigned to library / read group %s ...\n" "${read_file_name}" "${sample_flowcell}" "${sample_lane}" "${sample_name}" "${RG_LB}"
    printf "\nCreating human bam file: %s.sort.human.bam\n" "${bam_dir}/${read_file_name}"

    bwa mem \
        -t "${n_threads}" \
        -Y -K 100000000 \
        -R "@RG\tID:${RG_ID}\tSM:${RG_SM}\tPL:ILLUMINA\\tPU:${RG_PU}\\tLB:${RG_LB}" \
        "${ref_human}" \
        "${read_file_path}_R1_001.trim.fastq.gz" \
        "${read_file_path}_R2_001.trim.fastq.gz" |
    # sort and compress to bam
        samtools sort --threads "${n_threads}" \
            -o "${bam_dir}/${read_file_name}.sort.human.bam"

    # extract all unmapped pairs (both reads unmapped)
    # approach adapted from https://lh3.github.io/2021/07/06/remapping-an-aligned-bam
    # TODO: alternatively use bedtools' bamtofastq approach and save intermediate steps
    printf "\nMapping human filtered reads to %s genome for sample %s...\n" "${species}" "${read_file_name}"
    printf "\nCreating bam file: %s.sort.bam\n" "${bam_dir}/${read_file_name}"
    samtools view -b -f 12 "${bam_dir}/${read_file_name}.sort.human.bam" |
    # convert back to fastq
        samtools collate -Oun128 - |
        samtools fastq -OT RG,BC - |
    # map to plasmodium genome
    # -CH adds back original read group info
    # -p gathers paired reads from stream - https://github.com/samtools/samtools/issues/1306
        bwa mem \
            -t "${n_threads}" \
            -Y -K 100000000 \
            -CH <(samtools view -H "${bam_dir}/${read_file_name}.sort.human.bam" | grep ^@RG) \
            -p \
            "${ref}" \
            - |
    # sort and compress to bam
        samtools sort --threads "${n_threads}" \
            -o "${bam_dir}/${read_file_name}.sort.bam"
    printf "\n----------------\nFinished aligning reads in ${read_file_path} R1/R2."
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
    for i in "${1}_"*".sort.bam"; do
        # echo "looping over ${i}"
        arr+=( "--INPUT ${i}" )
    done;
    echo ${arr[@]}
}
export -f add_input_prefix

# ! Note that single quotes are required to avoid the command substitution around the
# ! add_input_prefix function call from being executed before it is passed to parallel

if [[ "${fastq_identifier}" == "name_first" ]]; then
    field=1
elif [[ "${fastq_identifier}" == "flowcell_first" ]]; then
    field=2
fi

for i in "${bam_dir}/"*.sort.bam; do
    echo "$(basename "${i}" | cut -d '_' -f ${field} )";
done \
    | sort -u \
    | parallel -j "${jobs}" --halt now,fail=1 \
        gatk --java-options -Xmx${mem}G \
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
# rm "${bam_dir}/"*.sort.bam "${bam_dir}/"*.sort.human.bam

# aggregate results with multiQC
multiqc --force "${output_dir}" --config "${multiqc_conf}" --outdir "${output_dir}/multiqc"

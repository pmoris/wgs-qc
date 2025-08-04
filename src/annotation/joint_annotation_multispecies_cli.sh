#!/usr/bin/env bash
# Description: Script to perform VCF annotation using snpEff jointly on provided species
# Author: Pieter Moris

# set bash strict mode - optionally add x to show commands
# set -euo pipefail
# disabled because of the many gotchas, see https://mywiki.wooledge.org/BashPitfalls?highlight=%28pipefail%29#set_-euo_pipefail
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
Usage: ${0##*/} [-h]  [-o OUTPUT DIRECTORY ] [-s pf pv]
    -h                                      display this help and exit
    -o | --output_dir OUTPUT DIRECTORY      File path to output directory; should already
                                            contain joint VCF files.
                                            (default = PROJECT_ROOT/results/)
    -s | --species_list                     Species to perform joint variant calling for
                                            Options are: pf, pv, pm poc, pow, pk
EOF
}

while :; do
    case ${1:-} in
        -h|-\?|--help)
            show_help    # Display a usage synopsis.
            exit
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

        -s|--species_list)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                species_list=$2
                shift
            else
                die 'ERROR: "--species_list" requires a non-empty option argument.'
            fi
            ;;
        --species_list=?*)
            species_list=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --species_list=)         # Handle the case of an empty --output_dir=
            die 'ERROR: "--species_list" requires a non-empty option argument.'
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

if [ -n "${SLURM_MEM_PER_NODE-}" ]; then
    mem=$((${SLURM_MEM_PER_NODE}/1000}))
elif [ -n "${SLURM_MEM_PER_CPU-}" ]; then
    mem=$((${SLURM_MEM_PER_CPU}/1000*${n_threads}))
else
    mem=8
fi

# define in and outputs
output_dir="$(realpath "${output_dir:-"${PROJECT_ROOT}/results/"}")"
vcf_dir="${output_dir}/gatk/"
ann_dir="${output_dir}/snpeff/"

# create output directories
mkdir -p "${ann_dir}"

# reference database names
ref_pf="PlasmoDB-68_Pfalciparum3D7"
ref_pv="PlasmoDB-68_PvivaxPAM"
ref_pm="PlasmoDB-68_PmalariaeUG01"
ref_poc="PlasmoDB-68_PovalecurtisiGH01"
ref_pow="PlasmoDB-68_PovalewallikeriPowCR01"
ref_pk="PlasmoDB-68_PknowlesiH"

# config files
multiqc_conf="${PROJECT_ROOT}/config/multiqc_config.yaml"
snpeff_db="${PROJECT_ROOT}/data/snpEff_database/"
snpeff_conf="${PROJECT_ROOT}/config/snpEff.config"

# check if vcf directory exist
if [ ! -d "${vcf_dir}" ]; then
    echo "VCF input directory (${vcf_dir}) does not exist."
    exit 1
fi

# check if snpeff directory exists
if ! [ -d "${snpeff_db}" ]; then
    echo "SnpEff database not found in expected location: (${snpeff_db})."
    exit 1
fi

# log run options
printf "
Joint annotation script | $(basename "$0")
==============================================

Output directory:           ${ann_dir}
VCF directory:              ${vcf_dir}
Reference:                  ${snpeff_db}
Threads:                    ${n_threads}
Memory:                     ${mem}
"

#######################
# Start of annotation #
#######################

# for species in pv pm pk; do
for species in ${species_list}; do

    printf "\nAnnotating combined and filtered VCF files for ${species} using ${mem}G of memory...\n"

    # create output directories
    mkdir -p "${ann_dir}/${species}"

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
    elif [[ "${species}" == "pk" ]]; then
        ref="${ref_pk}"
    else
        printf "\n Provided species ${species} not supported. Please use any of pf|pv|pm|pow|poc|pk."
        exit 1
    fi

    # # create snpeff database
    # snpEff build -c "${snpeff_conf}" -gff3 -v -noCheckProtein "${ref}"

    # run snpeff annotation on combined vcf file
    snpEff -Xmx${mem}g -c "${snpeff_conf}" -v "${ref}" \
        -s "${ann_dir}/${species}/snpEff_summary_filter_added.html" \
        -csvStats "${ann_dir}/${species}/snpEff_summary_filter_added.csv" \
        "${vcf_dir}/${species}/combined.filter_added.vcf.gz" > "${ann_dir}/${species}/combined.filter_added.ann.vcf"

    snpEff -Xmx${mem}g -c "${snpeff_conf}" -v "${ref}" \
        -s "${ann_dir}/${species}/snpEff_summary_filtered_filtered.html" \
        -csvStats "${ann_dir}/${species}/snpEff_summary_filtered.csv" \
        "${vcf_dir}/${species}/combined.filtered.vcf.gz" > "${ann_dir}/${species}/combined.filtered.ann.vcf"

    gatk VariantsToTable -V "${ann_dir}/${species}/combined.filter_added.ann.vcf" -F CHROM -F POS -F TYPE -GF GT -O "${ann_dir}/${species}/combined.filter_added.table"

    gatk VariantsToTable -V "${ann_dir}/${species}/combined.filtered.ann.vcf" -F CHROM -F POS -F TYPE -GF GT -O "${ann_dir}/${species}/combined.filtered.table"

    # aggregate results with multiQC
    printf "\nRunning MultiQC on ${output_dir}...\n"
    multiqc --force "${output_dir}" --config "${multiqc_conf}" --outdir "${output_dir}/multiqc"

done

printf "\n#######################\nEnd of joint annotation script\n#######################\n"

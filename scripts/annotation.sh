#!/usr/bin/env bash
# Description: Script to perform VCF annotation using snpEff
# Author: Pieter Moris

# TODO: check or remove single_species option

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
Usage: ${0##*/} [-h] [-s SAMPLESHEET.CSV ] [-o OUTPUT DIRECTORY ] [-p]
    -h                                      display this help and exit
    -s | --samplesheet SAMPLESHEET.CSV      File path to samplesheet with sample-species info
    -o | --output_dir OUTPUT DIRECTORY      File path to output directory; should already
                                            contain joint VCF files.
                                            (default = PROJECT_ROOT/results/)
    -p | --single_species                   Species override for all samples. Options are:
                                            pf, pv, pm poc, pow
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

        -p|--single_species)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                single_species=$2
                shift
            else
                die 'ERROR: "--single_species" requires a non-empty option argument.'
            fi
            ;;

        --single_species=?*)
            single_species=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --single_species=)         # Handle the case of an empty --output_dir=
            die 'ERROR: "--single_species" requires a non-empty option argument.'
            ;;

        --snpeff_config)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                snpeff_config=$2
                shift
            else
                die 'ERROR: "--snpeff_config" requires a non-empty option argument.'
            fi
            ;;
        --snpeff_config=?*)
            snpeff_config=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --snpeff_config=)         # Handle the case of an empty --samplesheet=
            die 'ERROR: "--samplesheet" requires a non-empty option argument.'
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
samplesheet="$(realpath "${samplesheet:-"${PROJECT_ROOT}/data/samplesheet.csv"}")"
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
# snpeff_conf="${PROJECT_ROOT}/config/snpEff.config"    # set this via cli option instead
# snpeff_db="${PROJECT_ROOT}/data/snpEff_database/"     # not used but could be set via -dataDir if needed (relative path from snpEff.config)

# check if samplesheet exist
if [ ! -f "${samplesheet}" ]; then
    die "Samplesheet (${samplesheet}) does not exist."
fi

# check if vcf directory exist
if [ ! -d "${vcf_dir}" ]; then
    die "VCF input directory (${vcf_dir}) does not exist."
fi

# check if snpeff config file exists
if ! [ -d "${snpeff_config}" ]; then
    die "SnpEff config file not found in expected location: (${snpeff_config})."
fi

# # check if snpeff directory exists
# if ! [ -d "${snpeff_db}" ]; then
#     die "SnpEff database not found in expected location: (${snpeff_db})."
# fi

# log run options
printf "
BWA MEM script | $(basename "$0")
==============================================

Output directory:           ${ann_dir}
VCF directory:              ${vcf_dir}
Samplesheet:                ${samplesheet}
Reference:                  ${snpeff_db}
Threads:                    ${n_threads}
Memory:                     ${mem}
"

#######################
# Start of annotation #
#######################

for species in $(tail -n+2 "${samplesheet}" | cut -f2 -d, | sort | uniq); do

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
    elif [[ -z "${single_species:-}" && "${single_species}" =~ ^(pf|pv|pm|pow|poc)$ ]]; then
        ref="${single_species}"
    fi
    if [ -z "${ref:-}" ]; then
        die "Unexpected species name ${species} found in samplesheet ${samplesheet} (or passed via --single_species). Exiting..."
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

    # # aggregate results with multiQC
    # multiqc --force "${output_dir}" --config "${multiqc_conf}" --outdir "${output_dir}/multiqc"

done


# alternative approach: loop through directories and find all species that way

# for combined_vcf in "${vcf_dir}/${species}/combined.filtered.vcf.gz"

# # unset species to make sure there are no left overs from previous iterations
# species=
# species=$(awk -v pat="${sample_id}" -F',' '$1 ~ pat { print $2; exit}' "${samplesheet}")
# if [ -z "${species:-}" ]; then
#     printf "\nCould not find sample ${bam} (search query = ${sample_id}) during species lookup in samplesheet ${samplesheet}. Exiting...\n"
#     exit 1
# fi

printf "\n#######################\nEnd of annotation script\n#######################\n"

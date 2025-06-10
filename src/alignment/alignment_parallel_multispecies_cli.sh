#!/usr/bin/env bash
# Description: Script to perform (trimmed) read alignment (competitive or filter-based) on Plasmodium and human host
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
Usage: ${0##*/} [-h] [-s SAMPLESHEET.CSV ] [-o OUTPUT DIRECTORY ]
                [-r1 READ 1 SUFFIX ] [-r2 READ 2 SUFFIX ] [-e READ FILE EXTENSION ]
                [-n <name_first/flowcell_first ]
    -h                                      display this help and exit
    -s | --samplesheet SAMPLESHEET.CSV      File path to samplesheet with sample-species info
    -o | --output_dir OUTPUT DIRECTORY      File path to output directory; should already
                                            contain trimmed reads directory named fastp
                                            (default = PROJECT_ROOT/results/)
    -r1 | --read_1_suffix R1_001            Suffix for read pair 1 (excluding file extension)
    -r2 | --read_2_suffix R2_001            Suffix for read pair 2 (excluding file extension)
    -e | --read_file_extension .fastq.gz    Read file extension
    -n | --fastq_identifier <string>        Specifies structure of fastq file name. Either
                                            "name_first", "flowcell_first", "novogene",
                                            "name_middle" or "SRA".
    -p | --single_species                   Species override for all samples. Options are:
                                            pf, pv, pm poc, pow
    -c | --competitive                      Enable competitive mapping mode instead of
                                            first filtering against the human genome.
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

        -c|--competitive)
            competitive="true"
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
if ! [[ "${fastq_identifier}" == "name_first" || "${fastq_identifier}" == "flowcell_first" || "${fastq_identifier}" == "novogene" || "${fastq_identifier}" == "name_middle" || "${fastq_identifier}" == "SRA" ]]; then
    printf "\nFastq identifier structure was not set correctly, please specify "name_first", "flowcell_first", "novogene", "name_middle" or "SRA".\n"
    exit 1
fi

# set alignment mode
if [[ -z "${competitive:-}" ]]; then
    alignment_mode="filter"
elif [[ "${competitive:-}" == "true" ]]; then
    alignment_mode="competitive"
else
    printf "\nAlignment mode not set correctly. --competitive is a standalone option. Omitting it uses filter-based mode.\n"
    exit 1
fi

# config files
multiqc_conf="${PROJECT_ROOT}/config/multiqc_config.yaml"

# reference genome files
ref_human="${PROJECT_ROOT}/data/ref/human/GRCh38.p14/GENCODE/47/GRCh38.primary_assembly.genome.fa.gz"
ref_pf="${PROJECT_ROOT}/data/ref/Pfalciparum/3D7/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"
ref_pv="${PROJECT_ROOT}/data/ref/Pvivax/PvPAM/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.fasta"
ref_pm="${PROJECT_ROOT}/data/ref/Pmalariae/UG01/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PmalariaeUG01_Genome.fasta"
ref_poc="${PROJECT_ROOT}/data/ref/Povale/curtisi/PocGH01/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PovalecurtisiGH01_Genome.fasta"
ref_pow="${PROJECT_ROOT}/data/ref/Povale/wallikeri/PowCR01/PlasmoDB-release-68/PlasmoDB-68_PovalewallikeriPowCR01_Genome.fasta"
ref_pk="${PROJECT_ROOT}/data/ref/Pknowlesi/H/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PknowlesiH_Genome.fasta"
ref_phix="${PROJECT_ROOT}/data/ref/PhiX/PhiX-NC_001422.1.fasta"

# combined reference genomes and bed files for competitive mapping

# # RefSeq
# ref_pf_combined="${PROJECT_ROOT}/data/ref/combined/human_pf/GRCh38.p14_3D7/RefSeq/GCF_000001405.40_GCF_000002765.6/concat-GCF_000001405.40_GRCh38.p14_genomic-GCF_000002765.6_GCA_000002765_genomic.fna.gz"
# ref_pv_combined="${PROJECT_ROOT}/data/ref/combined/human_pv/GRCh38.p14_PvPAM/RefSeq-GenBank/GCF_000001405.40_GCA_949152365.1/concat-GCF_000001405.40_GRCh38.p14_genomic-GCA_949152365.1_PVPAM_genomic.fna.gz"
# ref_pm_combined="${PROJECT_ROOT}/data/ref/combined/human_pm/GRCh38.p14_UG01/RefSeq/GCF_000001405.40_GCF_900090045.1/concat-GCF_000001405.40_GRCh38.p14_genomic-GCF_900090045.1_PmUG01_genomic.fna.gz"
# ref_poc_combined="${PROJECT_ROOT}/data/ref/combined/human_poc/GRCh38.p14_PocGH01/RefSeq-GenBank/GCF_000001405.40_GCA_900090035.2/concat-GCF_000001405.40_GRCh38.p14_genomic-GCA_900090035.2_PocGH01_genomic.fna.gz"
# ref_pow_combined="${PROJECT_ROOT}/data/ref/combined/human_pow/GRCh38.p14_PowCR01/RefSeq-Genbank/GCF_000001405.40_GCA_900090025.2/concat-GCF_000001405.40_GRCh38.p14_genomic-GCA_900090025.2_PowCR01_genomic.fna.gz"
# bed_pf="${PROJECT_ROOT}/data/ref/Pfalciparum/3D7/RefSeq/GCF_000002765.6/GCF_000002765.6_GCA_000002765_genomic.bed"
# bed_pv="${PROJECT_ROOT}/data/ref/Pvivax/PvPAM/GenBank/GCA_949152365.1_PVPAM/GCA_949152365.1_PVPAM_genomic.bed"
# bed_pm="${PROJECT_ROOT}/data/ref/Pmalariae/UG01/RefSeq/GCF_900090045.1/GCF_900090045.1_PmUG01_genomic.bed"
# bed_poc="${PROJECT_ROOT}/data/ref/Povale/curtisi/PocGH01/GenBank/GCA_900090035.2/GCA_900090035.2_PocGH01_genomic.bed"
# bed_pow="${PROJECT_ROOT}/data/ref/Povale/wallikeri/PowCR01/GenBank/GCA_900090025.2/GCA_900090025.2_PowCR01_genomic.bed"
# bed_pk="${PROJECT_ROOT}/data/ref/Pknowlesi/H/RefSeq/GCF_000006355.2/GCF_000006355.2_GCA_000006355.2_genomic.bed"

## Gencode-PlasmoDB
ref_pf_combined="${PROJECT_ROOT}/data/ref/combined/human_pf/GRCh38.p14_3D7/GENCODE-PlasmoDB/GENCODE-47-PlasmoDB-release-68/concat-PlasmoDB-68_Pfalciparum3D7_Genome-GRCh38.primary_assembly.genome.fna.gz"
ref_pv_combined="${PROJECT_ROOT}/data/ref/combined/human_pv/GRCh38.p14_PvPAM/GENCODE-PlasmoDB/GENCODE-47-PlasmoDB-release-68/concat-PlasmoDB-68_PvivaxPAM_Genome-GRCh38.primary_assembly.genome.fna.gz"
ref_pm_combined="${PROJECT_ROOT}/data/ref/combined/human_pm/GRCh38.p14_UG01/GENCODE-PlasmoDB/GENCODE-47-PlasmoDB-release-68/concat-PlasmoDB-68_PmalariaeUG01_Genome-GRCh38.primary_assembly.genome.fna.gz"
ref_poc_combined="${PROJECT_ROOT}/data/ref/combined/human_poc/GRCh38.p14_PocGH01/GENCODE-PlasmoDB/GENCODE-47-PlasmoDB-release-68/concat-PlasmoDB-68_PovalecurtisiGH01_Genome-GRCh38.primary_assembly.genome.fna.gz"
ref_pow_combined="${PROJECT_ROOT}/data/ref/combined/human_pow/GRCh38.p14_PowCR01/GENCODE-PlasmoDB/GENCODE-47-PlasmoDB-release-68/concat-PlasmoDB-68_PovalewallikeriPowCR01_Genome-GRCh38.primary_assembly.genome.fna.gz"
ref_pk_combined="${PROJECT_ROOT}/data/ref/combined/human_pk/GRCh38.p14_H/GENCODE-PlasmoDB/GENCODE-47-PlasmoDB-release-68/concat-PlasmoDB-68_PknowlesiH_Genome-GRCh38.primary_assembly.genome.fna.gz"
bed_pf="${PROJECT_ROOT}/data/ref/Pfalciparum/3D7/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_Pfalciparum3D7_Genome.bed"
bed_pv="${PROJECT_ROOT}/data/ref/Pvivax/PvPAM/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.bed"
bed_pm="${PROJECT_ROOT}/data/ref/Pmalariae/UG01/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PmalariaeUG01_Genome.bed"
bed_poc="${PROJECT_ROOT}/data/ref/Povale/curtisi/PocGH01/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PovalecurtisiGH01_Genome.be"
bed_pow="${PROJECT_ROOT}/data/ref/Povale/wallikeri/PowCR01/PlasmoDB-release-68/PlasmoDB-68_PovalewallikeriPowCR01_Genome.bed"
bed_pk="${PROJECT_ROOT}/data/ref/Pknowlesi/H/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PknowlesiH_Genome.bed"

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
for ref in ${ref_human} ${ref_pf} ${ref_pv} ${ref_poc} ${ref_pow} ${ref_pm} ${ref_pk}; do
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
Reference Pknowlesi:        ${ref_pk}
threads:                    ${n_threads}
Alignment mode:             ${alignment_mode}
"

####################
# Start of mapping #
####################

# create reference index if it does not yet exist
for ref in ${ref_human} ${ref_pf} ${ref_pv} ${ref_pm} ${ref_pow} ${ref_poc} ${ref_pk}; do
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

    # get filepath containing basename of each read pair without R1/R2 suffix or file extension
    read_file_path="${r1%${read_1_suffix}.trim.fastq.gz}"
    # read_file_basename=$(basename "${read_1}" "_R1_001.trim.fastq.gz")

    # remove filepath prefix
    read_file_basename="${read_file_path##*/}"

    # TODO: modify this to read from samplesheet.csv instead
    # retrieve lane and sample group

    if [[ "${fastq_identifier}" == "name_first" ]]; then
        # ANT5797_S262_L001_R1_001.trim.fastq.gz
        # ERR5740747_1.trim.fastq.gz
        sample_name="$(echo "${read_file_basename}" | cut -d '_' -f1)"
        # sample="${read_file_basename%%_*}"
        sample_lane="$(echo "${read_file_basename}" | grep -Po '_L\d{3}' | sed 's/_//g' || echo "L001")"    # fall back on L001 if missing, to stop set -u from stopping script
        # sample_lane="${read_file_basename##*_}"
        sample_flowcell="$(zcat "${r1}" | head -n 1 | cut -d ':' -f3)" || true

        sample_library="${sample_name}"
        # it is unclear whether or not the S### identifier refers to unique libraries or not
        # sample_group="${read_file_basename%%_L*}"
        # sample_library="${sample_group##*_}"

    elif [[ "${fastq_identifier}" == "flowcell_first" ]]; then
        # 22NY35LT3_106264-002-098_CCTCCTTT-CTTTCGCG_L007_R2.fastq.gz

        sample_name="$(echo "${read_file_basename}" | cut -d '_' -f2)"
        sample_lane="$(echo "${read_file_basename}" | grep -Po '_L\d{3}' | sed 's/_//g')"
        sample_flowcell="$(echo "${read_file_basename}" | cut -d '_' -f1)"
        sample_library="${sample_name}"

    elif [[ "${fastq_identifier}" == "novogene" ]]; then
        # ANT_5975_WB_EKDN250004460-1A_22YCG2LT3_L8_1.fq
        # ANT_5975_WB_EKDN250004460-1A_22YCG2LT3_L8_1.trim.fastq.gz

        # sample_name="$(echo "${read_file_basename}" | awk -F '_' '{ NF=NF-3; print }' OFS='_' )"  # number of columns needs to be checked carefully
        # echo "ANT_5975_WB_EKDN250004460-1A_22YCG2LT3_L8" | awk -F '_' '{ NF=NF-3; print }' OFS='_'
        # echo "ANT_6099_WB_EKDN250004469-1A_22M5WWLT4_L6" | rev | cut -f4- -d '_' | rev
        # echo "ANT_6099_EKDN250004469-1A_22M5WWLT4_L6" | awk -F '_' '{ NF=NF-3; print }' OFS='_'
        # echo "ANT_6099_EKDN250004469-1A_22M5WWLT4_L6" | rev | cut -f4- -d '_' | rev
        # echo "6099_EKDN250004469-1A_22M5WWLT4_L6" | awk -F '_' '{ NF=NF-3; print }' OFS='_'
        # echo "6099_EKDN250004469-1A_22M5WWLT4_L6" | rev | cut -f4- -d '_' | rev
        sample_name="$(echo "${read_file_basename}" | rev | cut -f4- -d '_' | rev )"
        sample_lane="$(echo "${read_file_basename}" | grep -Po '_L\d{1}' | sed 's/_//g')"
        sample_flowcell="$(zcat "${r1}" | head -n 1 | cut -d ':' -f3)" || true
        sample_library="$(echo "${read_file_basename}" | rev | cut -f3 -d '_' | rev )"

    elif [[ "${fastq_identifier}" == "name_middle" ]]; then
        # FCHTVJYCCXY_L2_WHRDMALtbmRABCAA-111_1.fq.gz
        sample_name="$(echo "${read_file_basename}" | cut -d '_' -f3)"
        sample_lane="$(echo "${read_file_basename}" | grep -Po '_L\d{1,3}_' | sed 's/_//g')"
        sample_flowcell="$(zcat "${r1}" | head -n 1 | cut -d ':' -f3)" || true
        sample_library="${sample_name}"

    elif [[ "${fastq_identifier}" == "SRA" ]]; then
        # SAMEA1527532_ERX151701_ERR175555_1.fastq.gz
        sample_name=$(basename "${read_file_basename}" | cut -d '_' -f1)
        sample_lane=$(basename "${read_file_basename}" | cut -d '_' -f3)
        sample_flowcell="$(zcat "${r1}" | head -n 1 | cut -d ':' -f3)" || true
        sample_library=$(basename "${read_file_basename}" | cut -d '_' -f2)

    fi

    # for i in "${bam_dir}/"*.sort.bam; do
    #     sample_name=$(basename ${i} | awk -F '_' '{ NF=NF-3; print }' OFS='_' ); echo $(basename $i); echo ${sample_name};
    # done
    # ANT_5975_WB_EKDN250004460-1A_22YCG2LT3_L8.sort.bam
    # ANT_5975_WB
    # ANT_6099_EKDN250004469-1A_22M5WWLT4_L6.sort.bam
    # ANT_6099

    # for r1 in "${trimmed_fastq_dir}/"*_1.trim.fastq.gz; do
    #     read_file_path="${r1%_1.trim.fastq.gz}";
    #     read_file_basename="${read_file_path##*/}";
    #     sample_name="$(echo "${read_file_basename}" | awk -F '_' '{ NF=NF-3; print }' OFS='_' )"; echo ${sample_name} - basename ${read_file_basename};
    # done
    # ANT_5975_WB - basename ANT_5975_WB_EKDN250004460-1A_22YCG2LT3_L8
    # ANT_6099 - basename ANT_6099_EKDN250004469-1A_22M5WWLT4_L6

    RG_ID="${sample_name}.${sample_flowcell}.${sample_lane}" #.barcode?
    RG_SM="${sample_name}"
    RG_PU="${sample_flowcell}.${sample_lane}"
    RG_LB="${sample_name}.${sample_library}"

    # unset species to make sure there are no left overs from previous iterations
    species=
    species=$(awk -v pat="${sample_name}" -F',' '$1 ~ pat { print $2; exit}' "${samplesheet}")
    if [ -z "${species:-}" ]; then
        printf "\nCould not find sample ${read_file_basename} (search query = ${sample_name}) during species lookup in samplesheet ${samplesheet}. Exiting...\n"
        exit 1
    fi

    # regular filter-based mapping
    # if [[ -z "${competitive:-}" ]]; then
    if [[ "${alignment_mode}" == "filter" ]]; then

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
        elif [[ "${species}" == "pk" ]]; then
            ref="${ref_pk}"
        elif [[ -z "${single_species:-}" && "${single_species}" =~ ^(pf|pv|pm|pow|poc|pk)$ ]]; then
            ref="${single_species}"
        fi
        if [ -z "${ref:-}" ]; then
            printf "\nCould not find correct reference based on species lookup in samplesheet ${samplesheet} (or wrong option passed for --single-species) for sample ${read_file_basename}. Exiting...\n"
            exit 1
        fi

        printf "\nFilter-based alignment mode was set: mapping to human reference genome ${ref_human} followed by mapping unaligned reads to parasite reference genome ${ref}.\n"

        printf "\nMapping raw reads to human reference for read file %s, flowcell %s, lane %s of sample %s, assigned to library / read group %s ...\n" "${read_file_basename}" "${sample_flowcell}" "${sample_lane}" "${sample_name}" "${RG_LB}"

        # map to human reference genome first to remove host reads
        if [[ ! -f "${bam_dir}/${read_file_basename}.sort.human.bam" ]]; then

            printf "\nCreating human bam file: %s.sort.human.bam\n" "${bam_dir}/${read_file_basename}"

            bwa mem \
                -t "${n_threads}" \
                -Y -K 100000000 \
                -R "@RG\tID:${RG_ID}\tSM:${RG_SM}\tPL:ILLUMINA\tPU:${RG_PU}\tLB:${RG_LB}" \
                "${ref_human}" \
                "${read_file_path}${read_1_suffix}.trim.fastq.gz" \
                "${read_file_path}${read_2_suffix}.trim.fastq.gz" |
            # sort and compress to bam
                samtools sort --threads "${n_threads}" \
                    -o "${bam_dir}/${read_file_basename}.sort.human.bam"
        else
            # skip human alignment if .sort.human.bam file is already present
            printf "\nHuman aligned BAM file found for ${read_file_basename}, skipping...\n"
        fi

        # extract all unmapped pairs (both reads unmapped)
        # approach adapted from https://lh3.github.io/2021/07/06/remapping-an-aligned-bam
        if [[ ! -f "${bam_dir}/${read_file_basename}.sort.bam" ]]; then
            # TODO: alternatively use bedtools' bamtofastq approach and save intermediate steps
            printf "\nMapping human filtered reads to %s genome for sample %s...\n" "${species}" "${read_file_basename}"
            printf "\nCreating parasite bam file: %s.sort.bam\n" "${bam_dir}/${read_file_basename}"

            samtools view -b -f 12 "${bam_dir}/${read_file_basename}.sort.human.bam" |
            # convert back to fastq
                samtools collate -Oun128 - |
                samtools fastq -OT RG,BC - |
            # map to plasmodium genome
            # -CH adds back original read group info
            # -p gathers paired reads from stream - https://github.com/samtools/samtools/issues/1306
                bwa mem \
                    -t "${n_threads}" \
                    -Y -K 100000000 \
                    -CH <(samtools view -H "${bam_dir}/${read_file_basename}.sort.human.bam" | grep ^@RG) \
                    -p \
                    "${ref}" \
                    - |
            # sort and compress to bam
                samtools sort --threads "${n_threads}" \
                    -o "${bam_dir}/${read_file_basename}.sort.bam"
        else
            # skip parasite alignment if .sort.bam file is already present
            printf "\nParasite aligned BAM files found for ${read_file_basename}, skipping...\n"
        fi

    # competitive mapping
    # else
    elif [[ "${alignment_mode}" == "competitive" ]]; then

        # unset ref to make sure there are no left overs from previous ref loop
        ref=
        if [[ "${species}" == "pf" ]]; then
            ref="${ref_pf_combined}"
            bed="${bed_pf}"
        elif [[ "${species}" == "pv" ]]; then
            ref="${ref_pv_combined}"
            bed="${bed_pv}"
        elif [[ "${species}" == "pm" ]]; then
            ref="${ref_pm_combined}"
            bed="${bed_pm}"
        elif [[ "${species}" == "pow" ]]; then
            ref="${ref_pow_combined}"
            bed="${bed_pow}"
        elif [[ "${species}" == "poc" ]]; then
            ref="${ref_poc_combined}"
            bed="${bed_poc}"
        elif [[ "${species}" == "pk" ]]; then
            ref="${ref_pk_combined}"
            bed="${bed_pk}"
        fi

        if ! [ -f "${bed}" ]; then
            index_files_found=0
            printf "\nCreating region-level bed file for ${ref}...\n"
            awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' "${ref}.fai" > "${bed}"
        fi

        # map reads to combined reference genome
        printf "\nCompetitive alignment mode was set: mapping to concatenated reference genome ${ref} followed by extracting ${species} parasite reads from bam file.\n"

        if [[ ! -f "${bam_dir}/${read_file_basename}.sort.combined.bam" ]]; then
            printf "\nCreating combined bam file: %s.sort.combined.bam" "${bam_dir}/${read_file_basename}"
            bwa mem \
                -t "${n_threads}" \
                -Y -K 100000000 \
                -R "@RG\tID:${RG_ID}\tSM:${RG_SM}\tPL:ILLUMINA\tPU:${RG_PU}\tLB:${RG_LB}" \
                "${ref}" \
                "${read_file_path}${read_1_suffix}.trim.fastq.gz" \
                "${read_file_path}${read_2_suffix}.trim.fastq.gz" |
            # sort and compress to bam
            samtools sort --threads "${n_threads}" \
                -o "${bam_dir}/${read_file_basename}.sort.combined.bam"
        else
            # skip combined alignment if .sort.human.bam file is already present
            printf "\nCombined human-parasite BAM file found for ${read_file_basename}, skipping...\n"
        fi

        # extract plasmodium reads
        if [[ ! -f "${bam_dir}/${read_file_basename}.sort.bam" ]]; then
            printf "\nExtracting parasite reads from combined bam file: %s.sort.combined.bam" "${bam_dir}/${read_file_basename}"
            printf "\nCreating parasite bam file: %s.sort.bam\n" "${bam_dir}/${read_file_basename}"
            samtools view -b -h -L "${bed}" "${bam_dir}/${read_file_basename}.sort.combined.bam" > "${bam_dir}/${read_file_basename}.sort.bam"
        else
            # skip extraction if parasite bam file is already present
            printf "\nParasite BAM file found for ${read_file_basename}, skipping...\n"
        fi

    # fallback option, should never be reached
    else
        printf "\nAlignment mode not set correctly. --competitive is a standalone option. Omitting it uses filter-based mode.\n"
        exit 1
    fi

    ################
    # fix readgroups
    # printf "Fixing RG readgroups for ${bam_dir}/${read_file_basename}.sort.bam:"
    # printf "@RG\tID:${RG_ID}\tSM:${RG_SM}\tPL:ILLUMINA\\tPU:${RG_PU}\\tLB:${RG_LB}"

    # samtools addreplacerg -r "@RG\tID:${RG_ID}\tSM:${RG_SM}\tPL:ILLUMINA\\tPU:${RG_PU}\\tLB:${RG_LB}" -m overwrite_all -o "${bam_dir}/${read_file_basename}.sort.bam.fixed" "${bam_dir}/${read_file_basename}.sort.bam"

    # printf "Fixing RG readgroups for ${bam_dir}/${read_file_basename}.sort.human.bam:"

    # samtools addreplacerg -r "@RG\tID:${RG_ID}\tSM:${RG_SM}\tPL:ILLUMINA\\tPU:${RG_PU}\\tLB:${RG_LB}" -m overwrite_all -o "${bam_dir}/${read_file_basename}.sort.human.bam.fixed" "${bam_dir}/${read_file_basename}.sort.human.bam"

    # # rename old files
    # mv "${bam_dir}/${read_file_basename}.sort.bam" "${bam_dir}/${read_file_basename}.sort.bam.wrongRG"
    # mv "${bam_dir}/${read_file_basename}.sort.human.bam" "${bam_dir}/${read_file_basename}.sort.human.bam.wrongRG"
    # # rename new files
    # mv "${bam_dir}/${read_file_basename}.sort.bam.fixed" "${bam_dir}/${read_file_basename}.sort.bam"
    # mv "${bam_dir}/${read_file_basename}.sort.human.bam.fixed" "${bam_dir}/${read_file_basename}.sort.human.bam"
    ################

    printf "\nFinished aligning reads in ${read_file_path} R1/R2.\n----------------\n"
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

function run_markduplicates() {
    # check if combined markdup file already exists
    if [[ -f "${bam_dir}/${1}.sort.markdup.bam" ]]; then
        # skip extraction if parasite bam file is already present
        printf "\nDuplicate-marked BAM file found for ${1}, skipping...\n"
        # return 0;
    else
        printf "\nCombining and marking duplicates for sample ${1} using input files: $(add_input_prefix ${1}).\n"
        gatk --java-options -Xmx${mem}G \
            MarkDuplicates \
                $(add_input_prefix "${1}") \
                --OUTPUT "${bam_dir}/${1}.sort.markdup.bam" \
                --METRICS_FILE "${bam_dir}/${1}.markdup.metrics" \
                --REMOVE_DUPLICATES false
    fi
}

function add_input_prefix() {
# NOTE: final _ or .suffix needs to be present after id, to ensure names are not extended like S1 -> S11

    # set correct file name pattern
    if [[ "${fastq_identifier}" == "name_first" ]]; then
        # 23060404_HTK3CDMXY_L001.sort.bam
        # ERR5740747.sort.bam
        pattern="${1}[_.]*sort.bam"
    elif [[ "${fastq_identifier}" == "flowcell_first" ]]; then
        # 22GTGTLT4_106264-001-113_TTGTTGCA-GACGTCGT_L008.sort.bam
        pattern="*_${1}_*.sort.bam"
    elif [[ "${fastq_identifier}" == "novogene" ]]; then
        # ANT_6745_EKDN250004467-1A_22M5WWLT4_L6.sort.bam
        pattern="${1}_*.sort.bam"
    elif [[ "${fastq_identifier}" == "name_middle" ]]; then
        # FCHTVJYCCXY_L2_WHRDMALtbmRABCAA-111.sort.bam
        pattern="*_${1}.sort.bam"
    elif [[ "${fastq_identifier}" == "SRA" ]]; then
        # # SAMEA1527532_ERX151701_ERR175555.sort.bam
        pattern="${1}_*.sort.bam"
    fi

    declare -a arr=()
    # NOTE: do not enclose pattern variable in quotes, as this will prevent glob pathname expansion
    for i in "${bam_dir}/"${pattern}; do
    # for i in "${bam_dir}/"*"${1}_"*".sort.bam"; do
        arr+=( "--INPUT ${i}" )
    done;
    echo ${arr[@]}
}

# export functions and variables so that they are accessible by parallel
export bam_dir
export fastq_identifier
export mem
export -f add_input_prefix
export -f run_markduplicates

for i in "${bam_dir}/"*.sort.bam; do
    # parse file name
    if [[ "${fastq_identifier}" == "name_first" ]]; then
        # 23060404_HTK3CDMXY_L001.sort.bam
        # ERR5740747.sort.bam -> note that .sort.bam suffix needs to be removed, because splitting on _ would retain it as it is part of the last element
        sample_name=$(basename "${i%.sort.bam}" | cut -d '_' -f1)
    elif [[ "${fastq_identifier}" == "flowcell_first" ]]; then
        # 22GTGTLT4_106264-001-113_TTGTTGCA-GACGTCGT_L008.sort.bam
        sample_name=$(basename "${i}" | cut -d '_' -f2)
    elif [[ "${fastq_identifier}" == "novogene" ]]; then
        # ANT_6745_EKDN250004467-1A_22M5WWLT4_L6.sort.bam
        sample_name=$(basename "${i}" | awk -F '_' '{ NF=NF-3; print }' OFS='_' )   # note that there is one fewer field now that _R1 and _R2 have been trimmed
    elif [[ "${fastq_identifier}" == "name_middle" ]]; then
        # FCHTVJYCCXY_L2_WHRDMALtbmRABCAA-111.sort.bam
        # note that .sort.bam suffix needs to be removed, because splitting on _ would retain it as it is part of the last element
        sample_name=$(basename "${i%.sort.bam}" | cut -d '_' -f3)
    elif [[ "${fastq_identifier}" == "SRA" ]]; then
        # # SAMEA1527532_ERX151701_ERR175555.sort.bam
        sample_name=$(basename "${i}" | cut -d '_' -f1)
    fi
    # echo sample name to pass it to parallel
    echo "${sample_name}";
done \
    | sort -u \
    | parallel -j "${jobs}" --halt now,fail=1 \
        run_markduplicates {}

# for bam in results-testset/bwa/*.sort.bam; do echo "${bam%%_L*}"; done | sort -u | while read -r line ; do array=(${line}*.sort.bam); echo ${array[@]}; done

# └─▶ for bam in results-testset/bwa/*.sort.bam; do echo "${bam%%_L*}"; done | sort -u | while read -r line ; do echo $(cat "sdfgsd" ${line}*.sort.bam); done

# └─▶ for i in results-testset/bwa/*.sort.bam; do echo "${i%%_L*}"; done | sort -u | parallel echo "{}*.sort.bam" "{}"

printf "\nCollecting alignment stats on bam files...\n"

# TODO: incorporate into parallel?
for bam in "${bam_dir}/"*.sort.markdup.bam; do
    if [[ -f "${bam}.bai" && -f "${bam}.stats" && -f "${bam}.flagstat" && -f "${bam}.idxstats" ]]; then continue; fi
    printf "\nCreating index and samtool stats for ${bam}...\n"
    samtools index --threads "${n_threads}" "${bam}"
    samtools stats --threads "${n_threads}" "${bam}" >"${bam}.stats"
    samtools flagstat --threads "${n_threads}" "${bam}" >"${bam}.flagstat"
    samtools idxstats --threads "${n_threads}" "${bam}" >"${bam}.idxstats"
done

if [[ "${alignment_mode}" == "filter" ]]; then
    alignment_mode_bam_suffix="human"
elif [[ "${alignment_mode}" == "competitive" ]]; then
    alignment_mode_bam_suffix="combined"
fi
for bam in "${bam_dir}/"*.sort.${alignment_mode_bam_suffix}.bam; do
    if [[ -f "${bam}.bai" && -f "${bam}.stats" && -f "${bam}.flagstat" && -f "${bam}.idxstats" ]]; then continue; fi
    printf "\nCreating index and samtool stats for ${bam}...\n"
    samtools index --threads "${n_threads}" "${bam}"
    samtools stats --threads "${n_threads}" "${bam}" >"${bam}.stats"
    samtools flagstat --threads "${n_threads}" "${bam}" >"${bam}.flagstat"
    samtools idxstats --threads "${n_threads}" "${bam}" >"${bam}.idxstats"
done

# clean up
# rm "${bam_dir}/"*.sort.bam "${bam_dir}/"*.sort.${alignment_mode_bam_suffix}.bam

# aggregate results with multiQC
printf "\nRunning MultiQC on ${output_dir}...\n"
multiqc --force "${output_dir}" --config "${multiqc_conf}" --outdir "${output_dir}/multiqc"

printf "\n#######################\nEnd of alignment script\n#######################\n"

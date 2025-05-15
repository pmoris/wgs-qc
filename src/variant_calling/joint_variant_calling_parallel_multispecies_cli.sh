#!/usr/bin/env bash
# Description: Script to perform joint variant calling on bam files jointly on provided species
# Author: Pieter Moris

# TODO: set GATK -verbosity to WARNING or ERROR instead of default INFO
# TODO: automatically clean up genomicsDB workspace and tmp after processing each species or interval (latter option would require parallel to call custom function with separate GATK and removal steps)

# set bash strict mode - optionally add x to show commands
# set -euo pipefail
# disabled because of the many gotchas, see https://mywiki.wooledge.org/BashPitfalls?highlight=%28pipefail%29#set_-euo_pipefail

# allow debug mode by running `TRACE=1 ./script.sh`
if [[ "${TRACE-0}" == "1" ]]; then set -x; fi

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
    -i | --input_dir INPUT DIRECTORY        File path to input directory with .g.vcf files
                                            (default = PROJECT_ROOT/results/)
    -o | --output_dir OUTPUT DIRECTORY      File path to output directory
                                            (default = PROJECT_ROOT/results/)
    -r1 | --read_1_suffix R1_001            Suffix for read pair 1 (excluding file extension)
    -r2 | --read_2_suffix R2_001            Suffix for read pair 2 (excluding file extension)
    -e | --read_file_extension .fastq.gz    Read file extension
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

        -i|--input_dir)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                input_dir=$2
                shift
            else
                die 'ERROR: "--input_dir" requires a non-empty option argument.'
            fi
            ;;
        --input_dir=?*)
            input_dir=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --input_dir=)         # Handle the case of an empty --input_dir=
            die 'ERROR: "--input_dir" requires a non-empty option argument.'
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

# define default in and outputs
input_dir="$(realpath "${input_dir:-"${PROJECT_ROOT}/results/"}")"
output_dir="$(realpath "${output_dir:-"${PROJECT_ROOT}/results/"}")"
vcf_dir="${output_dir}/gatk/"

# output directories will be created per species later on in the script
# mkdir -p "${vcf_dir}" "${vcf_dir}/haplotypecaller" "${vcf_dir}/genomicsdbimport" "${vcf_dir}/genotypegvcfs" "${vcf_dir}/variantfilter/snp" "${vcf_dir}/variantfilter/indel"

# config files
multiqc_conf="${PROJECT_ROOT}/config/multiqc_config.yaml"

# reference files
# ref_human="${PROJECT_ROOT}/data/ref/human/GRCh38.p14/GENCODE/47/GRCh38.primary_assembly.genome.fa.gz"
ref_pf="${PROJECT_ROOT}/data/ref/Pfalciparum/3D7/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"
ref_pv="${PROJECT_ROOT}/data/ref/Pvivax/PvPAM/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.fasta"
ref_pm="${PROJECT_ROOT}/data/ref/Pmalariae/UG01/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PmalariaeUG01_Genome.fasta"
ref_poc="${PROJECT_ROOT}/data/ref/Povale/curtisi/PocGH01/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PovalecurtisiGH01_Genome.fasta"
ref_pow="${PROJECT_ROOT}/data/ref/Povale/wallikeri/PowCR01/PlasmoDB-release-68/PlasmoDB-68_PovalewallikeriPowCR01_Genome.fasta"
ref_pk="${PROJECT_ROOT}/data/ref/Pknowlesi/H/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PknowlesiH_Genome.fasta"
# ref_phix="${PROJECT_ROOT}/data/ref/PhiX/PhiX-NC_001422.1.fasta"

# intervals="/data/antwerpen/grp/ap_itg_mu/public_data/reference_genomes/Pvivax/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.bed"
# intervals should be a bed file containing only the plasmodium regions. Using the bed file for the
# full combined reference genome creates empty files and requires more resources.
# bed file will be created automatically or searched for using the full file name and replacing .fasta(.gz) with .bed

# TODO change combined reference genome and intervals files like in alignment.sh

# check if reference fasta files exists
for ref in ${ref_pf} ${ref_pv} ${ref_poc} ${ref_pow} ${ref_pm} ${ref_pk}; do
    if ! [ -f "${ref}" ]; then
        printf "\nReference fasta file not found (${ref}).\n"
        exit 1
    fi
done

# log run options
printf "
Joint variant call script | $(basename "$0")
==============================================

Output directory:           ${vcf_dir}


Reference Pfalciparum:      ${ref_pf}
Reference Pvivax:           ${ref_pv}
Reference Pmalaria:         ${ref_pm}
Reference Povale wallikeri: ${ref_pow}
Reference Povale curtisi:   ${ref_poc}
Reference Pknowlesi:        ${ref_pk}
threads:                    ${n_threads}
"

############################
# Start of variant calling #
############################

for species in ${species_list}; do

    # create output directories
    mkdir -p "${vcf_dir}/${species}" \
        "${vcf_dir}/${species}/haplotypecaller" \
        "${vcf_dir}/${species}/genomicsdbimport" \
        "${vcf_dir}/${species}/genotypegvcfs" \
        "${vcf_dir}/${species}/variantfilter/snp" \
        "${vcf_dir}/${species}/variantfilter/indel"

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

    # Create reference fai and dict files if they do not yet exist
    index_files_found=1
    if ! [ -f "${ref}.fai" ]; then
        index_files_found=0
        printf "\nBuilding samtools faidx index for ${ref}...\n"
        samtools faidx "${ref}"
    fi
    # if ! [ -f "${ref%.fasta.gz}.dict" ]; then
    if ! [ -f "${ref%.fasta}.dict" ]; then
        # NOTE: GATK CreateSequenceDictionary only accepts `ref.fasta(.gz)` files and
        # outputs `ref.dict` by default (when not using -O).
        # Downstream tools automatically look for `ref.dict` (i.e., you cannot supply
        # your own dictionary file). Thus, when looking for the existence of dict file,
        # 1) bgz or fna files need to be renamed
        # 2) full suffix needs to be stripped for referring to the dict file
        index_files_found=0
        printf "\nBuilding GATK sequence dictionary for ${ref}...\n"
        gatk CreateSequenceDictionary -R "${ref}"
    fi
    if ! [ -f "${ref}.bed" ]; then
        index_files_found=0
        printf "\nCreating region-level bed file for ${ref}...\n"
        awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' "${ref}.fai" > "${ref}.bed"
    fi
    if [ "$index_files_found" -eq 1 ]; then
        printf "\nFound fai, dict and bed files for ${ref}, skipping indexing steps...\n"
    fi
    #  TODO: check which extension is used for reference
    # && [ -f "${ref%.fa.gz}.dict" ] && [ -f "${ref%.fa}.dict" ] && [ -f "${ref%.fasta}.dict" ] && [ -f "${ref%.fasta.gz}.dict" ]
    # create sequence dictionary automatically outputs file as basename .dict, so to check for presence the exact name of the file should be known
done

# parallellize by contig/region
jobs=$((${n_threads}/4))
if [ -n "${SLURM_MEM_PER_NODE-}" ]; then
    mem=$((${SLURM_MEM_PER_NODE}/1000/${jobs}))
elif [ -n "${SLURM_MEM_PER_CPU-}" ]; then
    mem=$((${SLURM_MEM_PER_CPU}/1000*4))
else
    mem=4
fi

# Combine gvfcs for each species, perform joint genotyping and filter variants
# for species in pv pm pk; do
# for species in pm pk; do
for species in ${species_list}; do

    printf "\n----------------\nCombining GVCFS, performing joint genotyping and filtering variants for all ${species} samples.\n"

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

    # skip if combined.filtered.vcf already exists
    if [[ -f "${vcf_dir}/${species}/combined.filtered.vcf.gz" ]]; then
        printf "\nCombined filtered VCF files found for ${species}, skipping...\n"
        continue
    fi

    # set bed file with intervals/regions for reference species - expected file name is the same as the .fasta ref, but with a .bed extension
    intervals=
    intervals="${ref%.fasta}.bed"

    # Create tmp and cache directories for genomicsdbimport
    # Note: tmp should already exist, workspace cache cannot exist yet ("${vcf_dir}/workspace")
    mkdir -p "${vcf_dir}/${species}/genomicsdbimport/tmp"

    # Create sample maps for genomicsdbimport
    cut -f1 "${intervals}" | \
    while read -r region; do
        # create sample map per region
        > "${vcf_dir}/${species}/genomicsdbimport/sample.${region}.map"

        # add region-specific vcf file for each sample
        # for gvcf in "${vcf_dir}/${species}/haplotypecaller/"*.${region}.g.vcf.gz; do
        find "${input_dir}" -name "*${region}.g.vcf.gz" | grep "/${species}/" | while read -r gvcf; do
            echo "$(basename ${gvcf} .${region}.g.vcf.gz)"$'\t'"${gvcf}" >> "${vcf_dir}/${species}/genomicsdbimport/sample.${region}.map"
        done
    done

    # Combine gvcf files per region
    printf "\nCombining GVCFs using GATK GenomicsDBImport, parallellized by region across ${jobs} jobs and assigning each ${mem}G of memory...\n"

    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" --halt now,fail=1 \
        gatk --java-options "-Xmx${mem}g" GenomicsDBImport \
            --genomicsdb-workspace-path "${vcf_dir}/${species}/genomicsdbimport/workspace-{}/" \
            --sample-name-map "${vcf_dir}/${species}/genomicsdbimport/sample.{}.map" \
            --tmp-dir "${vcf_dir}/${species}/genomicsdbimport/tmp" \
            --intervals "{}" \
            --overwrite-existing-genomicsdb-workspace \
            --batch-size 50 \
            --genomicsdb-shared-posixfs-optimizations true

    # Joint genotyping per region
    printf "\nPerforming joint genotyping using GATK GenotypeGVCFs, parallellized by region across ${jobs} jobs and assigning each ${mem}G of memory...\n"

    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" --halt now,fail=1 \
        gatk --java-options "-Xmx${mem}g" GenotypeGVCFs \
        -R "${ref}" \
        -V "gendb://${vcf_dir}/${species}/genomicsdbimport/workspace-{}" \
        -O "${vcf_dir}/${species}/genotypegvcfs/combined.{}.vcf.gz"

    # rm -r "${vcf_dir}/${species}/genomicsdbimport/workspace-"*
    # find results/ -name "workspace-*" -delete

    # filter variants - process snp and indels separately
    # See: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
    # https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
    # https://gatk.broadinstitute.org/hc/en-us/articles/360037499012-I-am-unable-to-use-VQSR-recalibration-to-filter-variants

    printf "\nFiltering variants per region (snp and indels separately), parallellized across ${jobs} jobs and assigning each ${mem}G of memory...\n"

    # snp
    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" --halt now,fail=1 \
        gatk --java-options "-Xmx${mem}g" SelectVariants \
            -V "${vcf_dir}/${species}/genotypegvcfs/combined.{}.vcf.gz" \
            -select-type SNP \
            -O "${vcf_dir}/${species}/variantfilter/snp/combined.{}.snp.vcf.gz"

    cut -f1 "${intervals}" | \
    parallel -q -j "${jobs}" --halt now,fail=1 \
        gatk --java-options "-Xmx${mem}g" VariantFiltration \
            -V "${vcf_dir}/${species}/variantfilter/snp/combined.{}.snp.vcf.gz" \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -filter "SOR > 3.0" --filter-name "StrandOddsRatio+3" \
            -O "${vcf_dir}/${species}/variantfilter/snp/combined.{}.snp.filter_added.vcf.gz"

    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" --halt now,fail=1 \
        gatk --java-options "-Xmx${mem}g" SelectVariants \
            -V "${vcf_dir}/${species}/variantfilter/snp/combined.{}.snp.filter_added.vcf.gz" \
            -R "${ref}" \
            --exclude-filtered true \
            -O "${vcf_dir}/${species}/variantfilter/snp/combined.{}.snp.filtered.vcf.gz"

    # indel
    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" --halt now,fail=1 \
        gatk --java-options "-Xmx${mem}g" SelectVariants \
            -V "${vcf_dir}/${species}/genotypegvcfs/combined.{}.vcf.gz" \
            -select-type INDEL \
            -O "${vcf_dir}/${species}/variantfilter/indel/combined.{}.indel.vcf.gz"

    cut -f1 "${intervals}" | \
    parallel -q -j "${jobs}" --halt now,fail=1 \
        gatk --java-options "-Xmx${mem}g" VariantFiltration \
            -V "${vcf_dir}/${species}/variantfilter/indel/combined.{}.indel.vcf.gz" \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "FS > 200.0" --filter-name "FS200" \
            -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
            -filter "SOR > 10.0" --filter-name "StrandOddsRatio+10" \
            -O "${vcf_dir}/${species}/variantfilter/indel/combined.{}.indel.filter_added.vcf.gz"

    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" --halt now,fail=1 \
        gatk --java-options "-Xmx${mem}g" SelectVariants \
            -V "${vcf_dir}/${species}/variantfilter/indel/combined.{}.indel.filter_added.vcf.gz" \
            -R "${ref}" \
            --exclude-filtered true \
            -O "${vcf_dir}/${species}/variantfilter/indel/combined.{}.indel.filtered.vcf.gz"

    # combine snp and indel vcf files for both filter_added and filtered vcf files
    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" --halt now,fail=1 \
        gatk --java-options "-Xmx${mem}g" SortVcf \
            -I "${vcf_dir}/${species}/variantfilter/indel/combined.{}.indel.filter_added.vcf.gz" \
            -I "${vcf_dir}/${species}/variantfilter/snp/combined.{}.snp.filter_added.vcf.gz" \
            -O "${vcf_dir}/${species}/variantfilter/combined.{}.filter_added.vcf.gz"

    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" --halt now,fail=1 \
        gatk --java-options "-Xmx${mem}g" SortVcf \
            -I "${vcf_dir}/${species}/variantfilter/indel/combined.{}.indel.filtered.vcf.gz" \
            -I "${vcf_dir}/${species}/variantfilter/snp/combined.{}.snp.filtered.vcf.gz" \
            -O "${vcf_dir}/${species}/variantfilter/combined.{}.filtered.vcf.gz"

    # combine regions for filtered and unfiltered vcf files

    printf "\nCombining filtered per-region VCF files (snp and indels separately), parallellized across ${jobs} jobs and assigning each ${mem}G of memory...\n"

    declare -a input_vcf_array=()
    # NOTE: contigs/intervals must be supplied in the genomic order!
    # for i in "${vcf_dir}/variantfilter/combined."*".filtered.vcf.gz"; do
    #     input_vcf_array+=( "--INPUT ${interval}" )
    # done;
    while read -r interval; do
        # echo "--INPUT \"${vcf_dir}/variantfilter/combined.${interval}.filtered.vcf.gz\""
        # NOTE: be careful with quotes! Quoting the entire array element passes on paths with double //
        # to gatk, resulting in java.nio.file.NoSuchFileException.
        # Omitting the outer quotes and only quoting the bash variable+string part does work.
        input_vcf_array+=( --INPUT "${vcf_dir}/${species}/variantfilter/combined.${interval}.filtered.vcf.gz" )
    done < <(cut -f1 "${intervals}")
    gatk GatherVcfs ${input_vcf_array[@]} -O "${vcf_dir}/${species}/combined.filtered.vcf.gz"

    declare -a input_vcf_array=()
    # for i in "${vcf_dir}/variantfilter/combined."*".filter_added.vcf.gz"; do
        # input_vcf_array+=( "--INPUT ${interval}" )
    # done;
    while read -r interval; do
        input_vcf_array+=( --INPUT "${vcf_dir}/${species}/variantfilter/combined.${interval}.filter_added.vcf.gz" )
    done < <(cut -f1 "${intervals}")
    gatk GatherVcfs ${input_vcf_array[@]} -O "${vcf_dir}/${species}/combined.filter_added.vcf.gz"

done

# aggregate results with multiQC
printf "\nRunning MultiQC on ${output_dir}...\n"
multiqc --force "${output_dir}" --config "${multiqc_conf}" --outdir "${output_dir}/multiqc"

printf "\n#######################\nEnd of joint variant calling script\n#######################\n"

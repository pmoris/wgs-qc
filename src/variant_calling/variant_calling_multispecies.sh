#!/usr/bin/env bash

##################################################
# Script to perform variant calling on bam files #
##################################################

# set bash strict mode
set -euo pipefail

# allow debug mode by running `TRACE=1 ./script.sh`
if [[ "${TRACE-0}" == "1" ]]; then set -x; fi

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
bam_dir="${output_dir}/bwa/"
vcf_dir="${output_dir}/gatk/"

multiqc_conf="${PROJECT_ROOT}/config/multiqc_config.yaml"

# create output directories
# mkdir -p "${vcf_dir}" "${vcf_dir}/haplotypecaller" "${vcf_dir}/genomicsdbimport" "${vcf_dir}/genotypegvcfs" "${vcf_dir}/variantfilter/snp" "${vcf_dir}/variantfilter/indel"

# reference files
# ref_human="${PROJECT_ROOT}/data/ref/human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
ref_pf="${PROJECT_ROOT}/data/ref/Pfalciparum/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"
ref_pv="${PROJECT_ROOT}/data/ref/Pvivax/PlasmoDB-68_PvivaxPAM_Genome.fasta"
ref_pm="${PROJECT_ROOT}/data/ref/Pmalariae/PlasmoDB-68_PmalariaeUG01_Genome.fasta"
ref_poc="${PROJECT_ROOT}/data/ref/Povale/curtisi/Poc221.fasta"
ref_pow="${PROJECT_ROOT}/data/ref/Povale/wallikeri/Pow222.fasta"
ref_phix="${PROJECT_ROOT}/data/ref/PhiX/PhiX-NC_001422.1.fasta"

# intervals="/data/antwerpen/grp/ap_itg_mu/public_data/reference_genomes/Pvivax/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.bed"
# intervals should be a bed file containing only the plasmodium regions. Using the bed file for the
# full combined reference genome creates empty files and requires more resources.

# check if bam directory exist
if [ ! -d "${bam_dir}" ]; then
    echo "BAM input directory (${bam_dir}) does not exist."
    exit 1
fi

# # check if reference fasta files exist
# if ! [ -f "${ref}" ]; then
#     echo "Reference fasta file not found (${ref})."
#     exit 1
# fi
# check if reference fasta files exists
for ref in ${ref_pf} ${ref_pv} ${ref_poc} ${ref_pow} ${ref_pm}; do
    if ! [ -f "${ref}" ]; then
        printf "\nReference fasta file not found (${ref}).\n"
        exit 1
    fi
done



# log run options
printf "
BWA MEM script | $(basename "$0")
==============================================

Output directory:           ${vcf_dir}
BAM directory:              ${bam_dir}
Reference Pv:               ${ref}
threads:                    ${n_threads}
"

############################
# Start of variant calling #
############################

for species in $(tail -n+2 "${PROJECT_ROOT}/data/samplesheet.csv" | cut -f2 -d, | sort | uniq); do

    # create output directories
    mkdir -p "${vcf_dir}/${species}" "${vcf_dir}/${species}/haplotypecaller" "${vcf_dir}/${species}/genomicsdbimport" "${vcf_dir}/${species}/genotypegvcfs" "${vcf_dir}/${species}/variantfilter/snp" "${vcf_dir}/${species}/variantfilter/indel"

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
    # intervals="${ref}.bed"
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

# Call variants per sample
for bam in "${bam_dir}"/*.sort.markdup.bam; do
    printf "\nCalling variants per region for ${bam}...\nParallellizing across ${jobs} jobs and assigning each ${mem}G of memory...\n"

    # get filepath containing basename of each read pair
    sample_path="${bam%.sort.markdup.bam}"
    # convert to basename of each read without the filepath prefix
    sample_name="${sample_path##*/}"

    # unset species to make sure there are no left overs from previous iterations
    species=
    species=$(awk -v pat="${sample_name}" -F',' '$1 ~ pat { print $2; exit}' "${PROJECT_ROOT}/data/samplesheet.csv")
    if [ -z "${species:-}" ]; then
        echo "Could not find sample during species lookup in samplesheet ${samplesheet}. Exiting..."
        exit 1
    fi

    # unset ref to make sure there are no left overs from previous ref loop
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

    intervals="${ref%.fasta}.bed"

    mkdir -p "${vcf_dir}/${species}/haplotypecaller"

    # sample_name=$(basename ${bam} .sort.markdup.bam)

    printf "\nCalling haplotypes for ${species} sample ${sample_name} using reference ${ref} and intervals ${intervals}\n"

    # cut -f1 "${ref}.bed" | \
    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" \
        gatk --java-options "-Xmx${mem}g" HaplotypeCaller \
            -R "${ref}" \
            -I "${bam}" \
            -O "${vcf_dir}/${species}/haplotypecaller/${sample_name}.{}.g.vcf.gz" \
            --native-pair-hmm-threads 4 \
            --intervals {} \
            -ERC GVCF
done

for species in $(tail -n+2 "${PROJECT_ROOT}/data/samplesheet.csv" | cut -f2 -d, | sort | uniq); do

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

    intervals="${ref%.fasta}.bed"

    # Create tmp and cache directories for genomicsdbimport
    # Note: tmp should already exist, workspace cache cannot exist yet
    mkdir -p "${vcf_dir}/${species}/genomicsdbimport/tmp" # "${vcf_dir}/workspace"

    # Create sample maps for genomicsdbimport
    cut -f1 "${intervals}" | \
    while read -r region; do
        # create sample map per region
        > "${vcf_dir}/${species}/genomicsdbimport/sample.$region.map"

        # add region-specific vcf file for each sample
        for gvcf in "${vcf_dir}/${species}/haplotypecaller/"*.${region}.g.vcf.gz; do
            echo "$(basename ${gvcf} .${region}.g.vcf.gz)"$'\t'"${gvcf}" >> "${vcf_dir}/${species}/genomicsdbimport/sample.$region.map"
        done
    done

    # Combine gvcf files per region
    printf "\nCombining gvcf files per region, parallellized across ${jobs} jobs and assigning each ${mem}G of memory...\n"

    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" \
        gatk --java-options "-Xmx${mem}g" GenomicsDBImport \
            --genomicsdb-workspace-path "${vcf_dir}/${species}/genomicsdbimport/workspace-{}/" \
            --sample-name-map "${vcf_dir}/${species}/genomicsdbimport/sample.{}.map" \
            --tmp-dir "${vcf_dir}/${species}/genomicsdbimport/tmp" \
            --intervals "{}" \
            --overwrite-existing-genomicsdb-workspace \
            --batch-size 50 \
            --genomicsdb-shared-posixfs-optimizations true

    # Joint genotyping per region
    printf "\nPerforming joint genotyping per region, parallellized across ${jobs} jobs and assigning each ${mem}G of memory...\n"

    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" \
        gatk --java-options "-Xmx${mem}g" GenotypeGVCFs \
        -R "${ref}" \
        -V "gendb://${vcf_dir}/${species}/genomicsdbimport/workspace-{}" \
        -O "${vcf_dir}/${species}/genotypegvcfs/combined.{}.vcf.gz"

    # filter variants - process snp and indels separately
    # See: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
    # https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
    # https://gatk.broadinstitute.org/hc/en-us/articles/360037499012-I-am-unable-to-use-VQSR-recalibration-to-filter-variants

    printf "\nFiltering variants per region (snp and indels separately), parallellized across ${jobs} jobs and assigning each ${mem}G of memory...\n"

    # snp
    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" \
        gatk --java-options "-Xmx${mem}g" SelectVariants \
            -V "${vcf_dir}/${species}/genotypegvcfs/combined.{}.vcf.gz" \
            -select-type SNP \
            -O "${vcf_dir}/${species}/variantfilter/snp/combined.{}.snp.vcf.gz"

    cut -f1 "${intervals}" | \
    parallel -q -j "${jobs}" \
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
    parallel -j "${jobs}" \
        gatk --java-options "-Xmx${mem}g" SelectVariants \
            -V "${vcf_dir}/${species}/variantfilter/snp/combined.{}.snp.filter_added.vcf.gz" \
            -R "${ref}" \
            --exclude-filtered true \
            -O "${vcf_dir}/${species}/variantfilter/snp/combined.{}.snp.filtered.vcf.gz"

    # indel
    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" \
        gatk --java-options "-Xmx${mem}g" SelectVariants \
            -V "${vcf_dir}/${species}/genotypegvcfs/combined.{}.vcf.gz" \
            -select-type INDEL \
            -O "${vcf_dir}/${species}/variantfilter/indel/combined.{}.indel.vcf.gz"

    cut -f1 "${intervals}" | \
    parallel -q -j "${jobs}" \
        gatk --java-options "-Xmx${mem}g" VariantFiltration \
            -V "${vcf_dir}/${species}/variantfilter/indel/combined.{}.indel.vcf.gz" \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "FS > 200.0" --filter-name "FS200" \
            -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
            -filter "SOR > 10.0" --filter-name "StrandOddsRatio+10" \
            -O "${vcf_dir}/${species}/variantfilter/indel/combined.{}.indel.filter_added.vcf.gz"

    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" \
        gatk --java-options "-Xmx${mem}g" SelectVariants \
            -V "${vcf_dir}/${species}/variantfilter/indel/combined.{}.indel.filter_added.vcf.gz" \
            -R "${ref}" \
            --exclude-filtered true \
            -O "${vcf_dir}/${species}/variantfilter/indel/combined.{}.indel.filtered.vcf.gz"

    # combine snp and indel vcf files for both filter_added and filtered vcf files
    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" \
        gatk --java-options "-Xmx${mem}g" SortVcf \
            -I "${vcf_dir}/${species}/variantfilter/indel/combined.{}.indel.filter_added.vcf.gz" \
            -I "${vcf_dir}/${species}/variantfilter/snp/combined.{}.snp.filter_added.vcf.gz" \
            -O "${vcf_dir}/${species}/variantfilter/combined.{}.filter_added.vcf.gz"

    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" \
        gatk --java-options "-Xmx${mem}g" SortVcf \
            -I "${vcf_dir}/${species}/variantfilter/indel/combined.{}.indel.filtered.vcf.gz" \
            -I "${vcf_dir}/${species}/variantfilter/snp/combined.{}.snp.filtered.vcf.gz" \
            -O "${vcf_dir}/${species}/variantfilter/combined.{}.filtered.vcf.gz"

    # combine regions for filtered and unfiltered vcf files
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

    # # aggregate results with multiQC
    # multiqc --force "${output_dir}" --config "${multiqc_conf}" --outdir "${output_dir}/multiqc"

done

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
mkdir -p "${vcf_dir}" "${vcf_dir}/haplotypecaller" "${vcf_dir}/genomicsdbimport" "${vcf_dir}/genotypegvcfs" "${vcf_dir}/variantfilter/snp" "${vcf_dir}/variantfilter/indel"

# reference files
# TODO: either make symlinks from the reference genome location to the
# data/ref/ directory in the project: e.g.,
# `ln -s /data/antwerpen/grp/ap_itg_mu/public_data/reference_genomes/Pvivax/PlasmoDB-release-68/* ./data/ref/`
# or use the full path to the reference genomes
# ref="${PROJECT_ROOT}/data/ref/GCF_000001405.40_GRCh38.p14_genomic-PlasmoDB-68_PvivaxPAM_Genome.fa.bgz"
# ref="${PROJECT_ROOT}/data/ref/GENCODE.release45-GRCh38.primary_assembly.genome-PlasmoDB-68_PvivaxPAM_Genome.fa.bgz"
# ref="/data/antwerpen/grp/ap_itg_mu/public_data/reference_genomes/combined/human_pv/GCF_000001405.40_GRCh38.p14_genomic-PlasmoDB-68_PvivaxPAM_Genome.fa.bgz"
ref="/data/antwerpen/grp/ap_itg_mu/public_data/reference_genomes/combined/human_pv/GENCODE.release45-GRCh38.primary_assembly.genome-PlasmoDB-68_PvivaxPAM_Genome.fa.gz"
# ref="${PROJECT_ROOT}/data/ref/Pvivax/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.fasta"
# ref="/data/antwerpen/grp/ap_itg_mu/public_data/reference_genomes/Pvivax/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.fasta"
intervals="/data/antwerpen/grp/ap_itg_mu/public_data/reference_genomes/Pvivax/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.bed"
# intervals should be a bed file containing only the plasmodium regions. Using the bed file for the
# full combined reference genome creates empty files and requires more resources.

# check if bam directory exist
if [ ! -d "${bam_dir}" ]; then
    echo "BAM input directory (${bam_dir}) does not exist."
    exit 1
fi

# check if reference fasta files exist
if ! [ -f "${ref}" ]; then
    echo "Reference fasta file not found (${ref})."
    exit 1
fi

# check if intervals bed file exist
if ! [ -f "${intervals}" ]; then
    echo "Bed file with region/contig intervals not found (${intervals})."
    exit 1
fi

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

# Create reference fai and dict files if they do not yet exist
index_files_found=1
if ! [ -f "${ref}.fai" ]; then
    index_files_found=0
    printf "\nBuilding samtools faidx index for ${ref}...\n"
    samtools faidx "${ref}"
fi
if ! [ -f "${ref}.dict" ]; then
    index_files_found=0
    printf "\nBuilding GATK sequence dictionary for ${ref}...\n"
    gatk CreateSequenceDictionary -R "${ref}" -O "${ref}.dict"
fi
if ! [ -f "${ref}.bed" ]; then
    index_files_found=0
    printf "\nCreating region-level bed file for ${ref}...\n"
    awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${ref}.fai > ${ref}.bed
fi
if [ "$index_files_found" -eq 1 ]; then
    printf "\nFound fai, dict and bed files for ${ref}, skipping indexing steps...\n"
fi
#  TODO: check which extension is used for reference
# && [ -f "${ref%.fa.gz}.dict" ] && [ -f "${ref%.fa}.dict" ] && [ -f "${ref%.fasta}.dict" ] && [ -f "${ref%.fasta.gz}.dict" ]
# create sequence dictionary automatically outputs file as basename .dict, so to check for presence the exact name of the file should be known

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
for bam in "${bam_dir}"/*.sort.plasmodium.markdup.bam; do
    printf "\nCalling variants per region for ${bam}...\nParallellizing across ${jobs} jobs and assigning each ${mem}G of memory...\n"

    sample_name=$(basename ${bam} .sort.plasmodium.markdup.bam)

    # cut -f1 "${ref}.bed" | \
    cut -f1 "${intervals}" | \
    parallel -j "${jobs}" \
    	gatk --java-options "-Xmx${mem}g" HaplotypeCaller \
            -R "${ref}" \
            -I "${bam}" \
            -O "${vcf_dir}/haplotypecaller/${sample_name}.{}.g.vcf.gz" \
            --native-pair-hmm-threads 4 \
            --intervals {} \
            -ERC GVCF
done

# Create tmp and cache directories for genomicsdbimport
# Note: tmp should already exist, workspace cache cannot exist yet
mkdir -p "${vcf_dir}/genomicsdbimport/tmp" # "${vcf_dir}/workspace"

# Create sample maps for genomicsdbimport
cut -f1 "${intervals}" | \
while read -r line; do
    # create sample map per region
    > "${vcf_dir}/genomicsdbimport/sample.$line.map"

    # add region-specific vcf file for each sample
    for gvcf in "${vcf_dir}/haplotypecaller/"*.${line}.g.vcf.gz; do
        echo "$(basename ${gvcf} .${line}.g.vcf.gz)"$'\t'"${gvcf}" >> "${vcf_dir}/genomicsdbimport/sample.$line.map"
    done
done

# Combine gvcf files per region
printf "\nCombining gvcf files per region, parallellized across ${jobs} jobs and assigning each ${mem}G of memory...\n"

cut -f1 "${intervals}" | \
parallel -j "${jobs}" \
    gatk --java-options "-Xmx${mem}g" GenomicsDBImport \
        --genomicsdb-workspace-path "${vcf_dir}/genomicsdbimport/workspace-{}/" \
        --sample-name-map "${vcf_dir}/genomicsdbimport/sample.{}.map" \
        --tmp-dir "${vcf_dir}/genomicsdbimport/tmp" \
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
    -V "gendb://${vcf_dir}/genomicsdbimport/workspace-{}" \
    -O "${vcf_dir}/genotypegvcfs/combined.{}.vcf.gz"

# filter variants - process snp and indels separately
# See: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
# https://gatk.broadinstitute.org/hc/en-us/articles/360037499012-I-am-unable-to-use-VQSR-recalibration-to-filter-variants

printf "\nFiltering variants per region (snp and indels separately), parallellized across ${jobs} jobs and assigning each ${mem}G of memory...\n"

# snp
cut -f1 "${intervals}" | \
parallel -j "${jobs}" \
    gatk --java-options "-Xmx${mem}g" SelectVariants \
        -V "${vcf_dir}/genotypegvcfs/combined.{}.vcf.gz" \
        -select-type SNP \
        -O "${vcf_dir}/variantfilter/snp/combined.{}.snp.vcf.gz"

cut -f1 "${intervals}" | \
parallel -q -j "${jobs}" \
    gatk --java-options "-Xmx${mem}g" VariantFiltration \
        -V "${vcf_dir}/variantfilter/snp/combined.{}.snp.vcf.gz" \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -filter "SOR > 3.0" --filter-name "StrandOddsRatio+3" \
        -O "${vcf_dir}/variantfilter/snp/combined.{}.snp.filter_added.vcf.gz"

cut -f1 "${intervals}" | \
parallel -j "${jobs}" \
    gatk --java-options "-Xmx${mem}g" SelectVariants \
        -V "${vcf_dir}/variantfilter/snp/combined.{}.snp.filter_added.vcf.gz" \
        -R "${ref}" \
        --exclude-filtered true \
        -O "${vcf_dir}/variantfilter/snp/combined.{}.snp.filtered.vcf.gz"

# indel
cut -f1 "${intervals}" | \
parallel -j "${jobs}" \
    gatk --java-options "-Xmx${mem}g" SelectVariants \
        -V "${vcf_dir}/genotypegvcfs/combined.{}.vcf.gz" \
        -select-type INDEL \
        -O "${vcf_dir}/variantfilter/indel/combined.{}.indel.vcf.gz"

cut -f1 "${intervals}" | \
parallel -q -j "${jobs}" \
    gatk --java-options "-Xmx${mem}g" VariantFiltration \
        -V "${vcf_dir}/variantfilter/indel/combined.{}.indel.vcf.gz" \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -filter "SOR > 10.0" --filter-name "StrandOddsRatio+10" \
        -O "${vcf_dir}/variantfilter/indel/combined.{}.indel.filter_added.vcf.gz"

cut -f1 "${intervals}" | \
parallel -j "${jobs}" \
    gatk --java-options "-Xmx${mem}g" SelectVariants \
        -V "${vcf_dir}/variantfilter/indel/combined.{}.indel.filter_added.vcf.gz" \
        -R "${ref}" \
        --exclude-filtered true \
        -O "${vcf_dir}/variantfilter/indel/combined.{}.indel.filtered.vcf.gz"

# combine snp and indel vcf files for both filter_added and filtered vcf files
cut -f1 "${intervals}" | \
parallel -j "${jobs}" \
    gatk --java-options "-Xmx${mem}g" SortVcf \
        -I "${vcf_dir}/variantfilter/indel/combined.{}.indel.filter_added.vcf.gz" \
        -I "${vcf_dir}/variantfilter/snp/combined.{}.snp.filter_added.vcf.gz" \
        -O "${vcf_dir}/variantfilter/combined.{}.filter_added.vcf.gz"

cut -f1 "${intervals}" | \
parallel -j "${jobs}" \
    gatk --java-options "-Xmx${mem}g" SortVcf \
        -I "${vcf_dir}/variantfilter/indel/combined.{}.indel.filtered.vcf.gz" \
        -I "${vcf_dir}/variantfilter/snp/combined.{}.snp.filtered.vcf.gz" \
        -O "${vcf_dir}/variantfilter/combined.{}.filtered.vcf.gz"

# combine regions for filtered and unfiltered vcf files
declare -a input_vcf_array=()
for i in "${vcf_dir}/variantfilter/combined."*".filtered.vcf.gz"; do
    input_vcf_array+=( "--INPUT ${i}" )
done;
gatk GatherVcfs ${input_vcf_array[@]} -O "${vcf_dir}/combined.filtered.vcf.gz"

declare -a input_vcf_array=()
for i in "${vcf_dir}/variantfilter/combined."*".filter_added.vcf.gz"; do
    input_vcf_array+=( "--INPUT ${i}" )
done;
gatk GatherVcfs ${input_vcf_array[@]} -O "${vcf_dir}/combined.filter_added.vcf.gz"

# # aggregate results with multiQC
# multiqc --force "${output_dir}" --config "${multiqc_conf}" --outdir "${output_dir}/multiqc"

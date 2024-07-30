#! /bin/env bash

#################################################################
# Script to perform quality control and trimming of fastq reads #
#################################################################

# set bash strict mode
set -euo pipefail

# allow debug mode by running `TRACE=1 ./qc.sh`
if [[ "${TRACE-0}" == "1" ]]; then set -x; fi

# get file path of script
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
# Alternatively, use absolute paths, but this makes the script less portable.

#####################
# Options and paths #
#####################

# set number of threads for downstream tools
n_threads="${SLURM_CPUS_PER_TASK:-8}"

# define in and outputs
# by defining everything here, there is no need to `cd` to directories first
.# fastq_dir="${SCRIPT_DIR}/data/fastq/"
output_dir="${SCRIPT_DIR}/results-reseq/"

fastq_screen_conf="${SCRIPT_DIR}/config/fastq-screen.conf"
multiqc_conf="${SCRIPT_DIR}/config/multiqc_config.yaml"

ref_human="/data/antwerpen/grp/aitg/arosanas/public_data/reference_genomes/human/Homo_sapiens.GRCh38.p14.GENCODE.release45/GRCh38.primary_assembly.genome.fa.gz"
ref_pf="/data/antwerpen/grp/aitg/arosanas/public_data/reference_genomes/Pfalciparum/PlasmoDB-release-68/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"
ref_pv="/data/antwerpen/grp/aitg/arosanas/public_data/reference_genomes/Pvivax/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.fasta"
ref_poc="/data/antwerpen/grp/aitg/arosanas/public_data/reference_genomes/Povale/PlasmoDB-release-68/PlasmoDB-68_PovalecurtisiGH01_Genome.fasta"
ref_pow="/data/antwerpen/grp/aitg/arosanas/public_data/reference_genomes/Povale/PlasmoDB-release-68/PlasmoDB-68_PovalewallikeriPowCR01_Genome.fasta"
ref_pm="/data/antwerpen/grp/aitg/arosanas/public_data/reference_genomes/Pmalariae/PlasmoDB-release-68/PlasmoDB-68_PmalariaeUG01_Genome.fasta"
#ref_Pk=""
ref_phix="/data/antwerpen/grp/aitg/arosanas/public_data/reference_genomes/PhiX/PhiX-NC_001422.1.fasta"

# define fastq read suffix
read_1_suffix="_R1_001.fastq.gz"
read_2_suffix="_R2_001.fastq.gz"
read_file_extension=".fastq.gz"

# create output directories
mkdir -p "${output_dir}/fastqc" \
    "${output_dir}/fastq-screen" \
    "${output_dir}/fastp" \
    "${output_dir}/multiqc"

# check if fastq directory exist
if [ ! -d "${fastq_dir}" ]; then
    echo "FASTQ input directory (${fastq_dir}) does not exist."
    exit 1
fi

# check if reference fasta files exists
for ref in ${ref_human} ${ref_pf} ${ref_pv} ${ref_poc} ${ref_pow} ${ref_pm} ${ref_phix}; do
    if ! [ -f "${ref}" ]; then
        echo "Reference fasta file not found (${ref})."
        exit 1
    fi
done

# log run options
printf "
FastQ QC script | $(basename "$0")
==============================================

Output directory:           ${output_dir}
FASTQ reads directory:      ${fastq_dir}
Reference human:            ${ref_human}
References:                 ${ref_pf} \t ${ref_pv} \t ${ref_poc} \t ${ref_pow} \t ${ref_pm} \t ${ref_phix}
threads:                    ${n_threads}
"

###############
# Start of QC #
###############

# create reference index if it does not yet exist
# required for fastq-screen
for ref in ${ref_human} ${ref_pf} ${ref_pv} ${ref_poc} ${ref_pow} ${ref_pm} ${ref_phix}; do
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

# run initial quality control
# note that --threads option only works when providing
# multiple files, i.e. glob instead of loop
# (alternatively use parallel
# find *.fq | parallel -j 10 "fastqc {} --outdir ..."
# or find exec )
echo "Running FastQC prior to trimming..."
fastqc \
    --threads "${n_threads}" \
    --outdir "${output_dir}/fastqc" \
    "${fastq_dir}/"*.fastq.gz

# run fastq-screen (threads option inherited by bwa/bowtie)
for read in "${fastq_dir}/"*.fastq.gz; do
    fastq_screen \
        --threads "${n_threads}" \
        --aligner bwa \
        --conf "${fastq_screen_conf}" \
        --outdir "${output_dir}/fastq-screen" \
        "${read}"
done

# trim reads (quality and adapters)
# fastp can use up to 16 threads, but efficiency is higher around 2-4 (https://hpc.nih.gov/training/gatk_tutorial/preproc.html#preproc-trim)
fastp_threads=$(( ${n_threads} > 16 ? 16 : ${n_threads} ))

for read_1 in "${fastq_dir}/"*"${read_1_suffix}"; do
    sample_name=$(basename "${read_1}" "${read_1_suffix}")
    read_2="${fastq_dir}/${sample_name}${read_2_suffix}"

    echo "Running fastp for ${sample_name}..."

    fastp \
        --in1 "${read_1}" \
        --in2 "${read_2}" \
        --out1 "${output_dir}/fastp/${sample_name}${read_1_suffix%${read_file_extension}}.fastp.fastq.gz" \
        --out2 "${output_dir}/fastp/${sample_name}${read_2_suffix%${read_file_extension}}.fastp.fastq.gz" \
        --json "${output_dir}/fastp/${sample_name}.fastp.json" \
        --html "${output_dir}/fastp/${sample_name}.fastp.html" \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_tail \
        --qualified_quality_phred 20 \
        --length_required 15 \
        --thread "${fastp_threads}"
    # Additional options:
    # --adapter_fasta                  specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file (string [=])
    # --unpaired1                      for PE input, if read1 passed QC but read2 not, it will be written to unpaired1. Default is to discard it. (string [=])
    # --unpaired2                      for PE input, if read2 passed QC but read1 not, it will be written to unpaired2. If --unpaired2 is same as --unpaired1 (default mode), both unpaired reads will be written to this same file. (string [=])
    # --failed_out                     specify the file to store reads that cannot pass the filters. (string [=])
    # --trim_poly_x                    enable polyX trimming in 3' ends.
    # --merge                          for paired-end input, merge each pair of reads into a single read if they are overlapped. The merged reads will be written to the file given by --merged_out, the unmerged reads will be written to the files specified by --out1 and --out2. The merging mode is disabled by default.
done

# TODO parallel option
# module load parallel
# find "${fastq_dir}" -name "*fastq.gz" |
    # parallel -j $((${n_threads} / 2)) \
        # fastp \
        # --in1 {} \
        # --in2 {=s/${read_1_suffix}/${read_2_suffix}} \
        # --out1 "${output_dir}/fastp/${sample_name}${read_1_suffix%${read_file_extension}}.fastp.fastq.gz" \
        # --out2 "${output_dir}/fastp/${sample_name}${read_2_suffix%${read_file_extension}}.fastp.fastq.gz" \
        # --json "${output_dir}/fastp/${sample_name}.fastp.json" \
        # --html "${output_dir}/fastp/${sample_name}.fastp.html" \
        # --detect_adapter_for_pe \
        # --cut_front \
        # --cut_tail \
        # --qualified_quality_phred 20 \
        # --length_required 15 \
        # --thread "${fastp_threads}"


# re-run qc after trimming
echo "Re-running FastQC after trimming..."
fastqc \
    --threads "${n_threads}" \
    --outdir "${output_dir}/fastqc/" \
    "${output_dir}/fastp/"*".fastp.fastq.gz"
# "${output_dir}/fastp/${sample_name}${read_1_suffix}.fastp.fastq.gz" \
# "${output_dir}/fastp/${sample_name}${read_2_suffix}.fastp.fastq.gz"

# aggregate results with multiQC
multiqc "${output_dir}" --config "${multiqc_conf}" --outdir "${output_dir}/multiqc"

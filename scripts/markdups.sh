#!/usr/bin/env bash

# set bash strict mode
set -euo pipefail

# allow debug mode by running `TRACE=1 ./script.sh` - equivalent to `set -x`
if [[ "${TRACE-0}" == "1" ]]; then set -o xtrace; fi

# get file path of script
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

#####################
# Options and paths #
#####################

# set number of threads for downstream tools
n_threads="${SLURM_CPUS_PER_TASK:-8}"

bam_dir="${SCRIPT_DIR}/../results/bwa/"
output_dir="${SCRIPT_DIR}/../results/markdups/"
mkdir -p "${output_dir}"

# picard mark duplicates
jobs=$((${n_threads}/2))

parallel -j "${jobs}" \
    gatk --java-options -Xmx$((8))G \
        MarkDuplicates \
        --INPUT "{}" \
        --OUTPUT "${output_dir}/{/.}.markdup.bam" \
        --METRICS_FILE {.}.markdup.metrics \
        --REMOVE_DUPLICATES false \
    ::: "${bam_dir}"/*.sort.bam


# "-Xmx$((8))G -XX:-UsePerfData"
# gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \

# # picard MarkDuplicates \
#     -I "${sample_name}.bam" \
#     -O "${output_dir}/${sample_name}.dups.bam" \
#     -M "${output_dir}/${sample_name}.markeddups_metrics.txt"
# '-REMOVE_DUPLICATES false -VALIDATION_STRINGENCY LENIENT'
# 
for bam in "${output_dir}"/*.sort.markdup.bam; do
    samtools stats --threads "${n_threads}" "${bam}" >"${bam}.stats"
    samtools flagstat --threads "${n_threads}" "${bam}" >"${bam}.flagstat"
    samtools idxstats --threads "${n_threads}" "${bam}" >"${bam}.idxstats"
done

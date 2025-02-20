#!/usr/bin/env bash

#SBATCH --job-name=wgs-variant-calling  # create a short name for your job
#SBATCH --nodes=1                   # node count
#SBATCH --ntasks=1                  # total number of tasks across all nodes
#SBATCH --cpus-per-task=28          # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=3750M         # memory per cpu-core (4G is default)
#SBATCH --time=24:00:00             # total run time limit (HH:MM:SS)
#SBATCH --mail-type=BEGIN,END,FAIL  # send email when job begins, ends or fails
#SBATCH --mail-user=pmoris@itg.be
#SBATCH -o %j-%x-stdout.out
#SBATCH -e %j-%x-stderr.out

printf "Script execution started at $(date).\n"

# define required software
# source /data/antwerpen/203/vsc20380/miniforge3/etc/profile.d/conda.sh
# conda activate wgs
export PATH="/scratch/antwerpen/203/vsc20380/containers/wgs/wgs-container/bin:$PATH"

module load parallel

# Get file path of script and set project root.
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR=$(dirname "$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')")
else
    SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
fi

echo "Script dir = ${SCRIPT_DIR}"
PROJECT_ROOT=$(realpath "${SCRIPT_DIR}/../")
echo "Project root = ${PROJECT_ROOT}"

# qc
bash "${PROJECT_ROOT}/scripts/qc.sh" \
    --input_dir "${PROJECT_ROOT}/data/fastq/test" \
    --output_dir "${PROJECT_ROOT}/results/test" \
    -r1 "_R1_001" -r2 "_R2_001" -e ".fastq.gz"

# bash "${PROJECT_ROOT}/scripts/qc.sh" \
#     --input_dir "${PROJECT_ROOT}/data/fastq/batch_001" \
#     --output_dir "${PROJECT_ROOT}/results/batch_001" \
#     -r1 "_R1_001" -r2 "_R2_001" -e ".fastq.gz"

# bash "${PROJECT_ROOT}/scripts/qc.sh" \
#     --input_dir "${PROJECT_ROOT}/data/fastq/batch_002" \
#     --output_dir "${PROJECT_ROOT}/results/batch_002" \
#     -r1 "_R1" -r2 "_R2" -e ".fastq.gz"

# bash "${PROJECT_ROOT}/scripts/qc.sh" \
#     --input_dir "${PROJECT_ROOT}/data/fastq/batch_003" \
#     --output_dir "${PROJECT_ROOT}/results/batch_003" \
#     -r1 "_R1" -r2 "_R2" -e ".fastq.gz"

# alignment
bash "${PROJECT_ROOT}/scripts/alignment.sh" \
    --samplesheet "${PROJECT_ROOT}/data/fastq/test/samplesheet.csv" \
    --output_dir "${PROJECT_ROOT}/results/test" \
    -r1 "_R1_001" -r2 "_R2_001" -e ".fastq.gz" \
    --fastq_identifier name_first

# bash "${PROJECT_ROOT}/scripts/alignment.sh" \
#     --samplesheet "${PROJECT_ROOT}/data/fastq/batch_001/samplesheet.csv" \
#     --output_dir "${PROJECT_ROOT}/results/batch_001" \
#     -r1 "_R1_001" -r2 "_R2_001" -e ".fastq.gz" \
#     --fastq_identifier name_first

# bash "${PROJECT_ROOT}/scripts/alignment.sh" \
#     --samplesheet "${PROJECT_ROOT}/data/fastq/batch_002/samplesheet.csv" \
#     --output_dir "${PROJECT_ROOT}/results/batch_002" \
#     -r1 "_R1" -r2 "_R2" -e ".fastq.gz" \
#     --fastq_identifier flowcell_first

# bash "${PROJECT_ROOT}/scripts/alignment.sh" \
#     --samplesheet "${PROJECT_ROOT}/data/fastq/batch_003/samplesheet.csv" \
#     --output_dir "${PROJECT_ROOT}/results/batch_003" \
#     -r1 "_R1" -r2 "_R2" -e ".fastq.gz" \
#     --fastq_identifier flowcell_first

# variant calling
bash "${PROJECT_ROOT}/scripts/variant-calling.sh" \
    --samplesheet "${PROJECT_ROOT}/data/fastq/test/samplesheet.csv" \
    --output_dir "${PROJECT_ROOT}/results/test" \
    --fastq_identifier name_first

# bash "${PROJECT_ROOT}/scripts/variant-calling.sh" \
#     --samplesheet "${PROJECT_ROOT}/data/fastq/batch_001/samplesheet.csv" \
#     --output_dir "${PROJECT_ROOT}/results/batch_001" \
#     --fastq_identifier name_first

# bash "${PROJECT_ROOT}/scripts/variant-calling.sh" \
#     --samplesheet "${PROJECT_ROOT}/data/fastq/batch_002/samplesheet.csv" \
#     --output_dir "${PROJECT_ROOT}/results/batch_002" \
#     --fastq_identifier flowcell_first


# bash "${PROJECT_ROOT}/scripts/variant-calling.sh" \
#     --samplesheet "${PROJECT_ROOT}/data/fastq/batch_003/samplesheet.csv" \
#     --output_dir "${PROJECT_ROOT}/results/batch_003" \
#     --fastq_identifier flowcell_first

# bash "${PROJECT_ROOT}/scripts/annotation.sh"

printf "Script finished at $(date).\n"

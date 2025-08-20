#!/usr/bin/env bash

#SBATCH --job-name=wgs-combined-genotyping   # create a short name for your job
#SBATCH --nodes=1                            # node count
#SBATCH --ntasks=1                           # total number of tasks across all nodes
#SBATCH --cpus-per-task=28                   # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=24:00:00                      # total run time limit (HH:MM:SS)
#SBATCH --mail-type=BEGIN,END,FAIL           # send email when job begins, ends or fails
#SBATCH --mail-user=pmoris@itg.be
#SBATCH -o %j-%x-stdout.out
#SBATCH -e %j-%x-stderr.out

printf "Script execution started at $(date).\n"

# define required software
# source /data/antwerpen/203/vsc20380/miniforge3/etc/profile.d/conda.sh
# conda activate wgs
export PATH="/scratch/antwerpen/203/vsc20380/containers/wgs/wgs-container-lock/bin:$PATH"

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

# joint genotyping: genomicsdbimport combining + filtering
bash "${PROJECT_ROOT}/scripts/variant-calling-combine.sh" \
    --input_dir "${PROJECT_ROOT}/input_gvcf/" \
    --output_dir "${PROJECT_ROOT}/results-combined/" \
    -s "pf pv pm po"

# annotation
bash "${PROJECT_ROOT}/scripts/annotation-combine.sh" \
    --output_dir "${PROJECT_ROOT}/results-combined/" \
    -s "pf pv pm po"

printf "Script finished at $(date).\n"

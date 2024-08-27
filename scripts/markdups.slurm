#!/bin/bash -l

#SBATCH --job-name=summit_markdups  # create a short name for your job
#SBATCH --nodes=1                   # node count
#SBATCH --ntasks=1                  # total number of tasks across all nodes
#SBATCH --cpus-per-task=28          # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G            # memory per cpu-core (4G is default)
#SBATCH --time=12:00:00             # total run time limit (HH:MM:SS)
#SBATCH --mail-type=BEGIN,END,FAIL  # send email when job begins, ends or fails
#SBATCH --mail-user=pmoris@itg.be
#SBATCH -o slurm-stdout.%j.%x.out
#SBATCH -e slurm-stderr.%j.%x.out

printf "Script execution started at $(date).\n"

module load parallel

source /data/antwerpen/203/vsc20380/miniforge3/etc/profile.d/conda.sh
conda activate wbc_depletion

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

bash "${SCRIPT_DIR}/markdups.sh"

printf "Script finished at $(date).\n"

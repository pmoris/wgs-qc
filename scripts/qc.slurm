#!/bin/bash

#SBATCH --job-name=summit_qc     # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=24       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=06:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=pmoris@itg.be

printf "Script execution started at $(date).\n"

# define required software
# source /data/antwerpen/203/vsc20380/miniforge3/bin/activate wgs

source /data/antwerpen/203/vsc20380/miniforge3/etc/profile.d/conda.sh
conda activate wgs

bash qc.sh

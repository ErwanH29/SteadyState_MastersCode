#!/bin/bash
#SBATCH --job-name=SteadyStateBH
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user="ehochart@gmail.com"
#SBATCH --mail-type="ALL"
#SBATCH --time=17:00:00
#SBATCH --partition=cpu-medium 
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=20

# load modules (assuming you start from the default environment)
# we explicitly call the modules to improve reproducibility
# in case the default settings change
module load SciPy-bundle/2019.10-foss-2019b-Python-3.7.4
module load Python/3.7.2-GCCcore-8.2.0
module load AMUSE/13.1.0-foss-2018a-Python-3.6.4

echo "[$SHELL] #### Starting Python test"
echo "[$SHELL] ## This is $SLURM_JOB_USER and this job has the ID $SLURM_JOB_ID"
# get the current working directory
export CWD=$(pwd)
echo "[$SHELL] ## current working directory: "$CWD

# Run the file
echo "[$SHELL] ## Run script"
python3 interface.py
echo "[$SHELL] ## Script finished"

echo "[$SHELL] #### Finished Python test. Have a nice day"

#!/bin/bash
#SBATCH --job-name=SteadyStateBH
#SBATCH --export=NONE
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
# module load pandas/0.24.0-foss-2018a-Python-3.6.4 
# module load pandas/0.24.0-foss-2019a-Python-3.7.2
module load AMUSE/2022.3-Miniconda3-4.11.0
# module load dotNET-Core-Runtime/5.0.4-GCCcore-9.3.0
# module load FFTW/3.3.8-gompic-2020b
# module load OpenMPI/4.0.5-gcccuda-2020b
# module load ScaLAPACK/2.1.0-gompic-2020b

echo "[$SHELL] #### Starting Script"
echo "[$SHELL] ## This is $SLURM_JOB_USER and this job has the ID $SLURM_JOB_ID"
# get the current working directory
export CWD=$(pwd)
echo "[$SHELL] ## current working directory: "$CWD

# Run the file
echo "[$SHELL] ## Run script"
python3 interface.py
echo "[$SHELL] ## Script finished"

echo "[$SHELL] #### Finished Python test. Have a nice day"

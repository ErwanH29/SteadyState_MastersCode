#!/bin/bash
#SBATCH --job-name=Testing_SteadyStateBH
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
<<<<<<< HEAD:Code/interface.slurm
#SBATCH --time=00:10:00
#SBATCH --partition=testing
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
=======
#SBATCH --mail-user="ehochart@gmail.com"
#SBATCH --mail-type="ALL"
#SBATCH --time=50:00:00
#SBATCH --partition=cpu-long 
#SBATCH --nodes=1
#SBATCH --ntasks=18
>>>>>>> 5e63cf24778c45ae8356dffcab987b50b519e175:ALICE_Comp_Scripts/interface.slurm

# load modules (assuming you start from the default environment)
# module load pandas/0.24.0-foss-2018a-Python-3.6.4 
<<<<<<< HEAD:Code/interface.slurm
# module load pandas/0.24.0-foss-2019a-Python-3.
#module load AMUSE/2022.3-Miniconda3-4.11.0
module load AMUSE/2022.7-Miniconda3-4.11.0
# source $HOME/python_test_venv/bin/activate
=======
# module load pandas/0.24.0-foss-2019a-Python-3.7.
module load AMUSE/2022.7-Miniconda3-4.11.0

# module load dotNET-Core-Runtime/5.0.4-GCCcore-9.3.0
# module load FFTW/3.3.8-gompic-2020b
# module load OpenMPI/4.0.5-gcccuda-2020b
# module load ScaLAPACK/2.1.0-gompic-2020b
>>>>>>> 5e63cf24778c45ae8356dffcab987b50b519e175:ALICE_Comp_Scripts/interface.slurm

export CWD=$(pwd)
export DATE=$(date)
echo "[$SHELL] #### Starting Script"
echo "[$SHELL] ## This is $SLURM_JOB_USER and this job has the ID $SLURM_JOB_ID"
# get the current working directory
echo "[$SHELL] ## current working directory: "$CWD

# Run the file
echo "[$SHELL] ## Run script"
<<<<<<< HEAD:Code/interface.slurm

python interface.py

echo "[$SHELL] ## Script finished"
echo "[$SHELL] ## Job done "$DATE
echo "[$SHELL] ## Used $SLURM_NTASKS cores"
echo "[$SHELL] ## Used $SLURM_CPUS_ON_NODE processors/CPUS"
echo "[$SHELL] ## Used $SLURM_CPUS_PER_TASK processors/CPUS per task"
=======
echo "[$SHELL] ## Job start "$DATE
python3 interface.py
echo "[$SHELL] ## Script finished"
echo "[$SHELL] ## Job done "$DATE
>>>>>>> 5e63cf24778c45ae8356dffcab987b50b519e175:ALICE_Comp_Scripts/interface.slurm

echo "[$SHELL] #### Finished Python test. Have a nice day"


#!/bin/bash
#SBATCH --partition=kamiak # Partition/Queue to use
#SBATCH --job-name=BEAST_Hanta # Job name
#SBATCH --output=BEAST_Hanta_S.out # Output file (stdout)
#SBATCH --error=BEAST_Hanta_S.err # Error file (stderr)
#SBATCH --mail-type=ALL # Email notification: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ricardo.rivero@wsu.edu # Email address for notifications
#SBATCH --time=7-00:00:00 # Wall clock time limit Days-HH:MM:SS
#SBATCH --nodes=1 # Number of nodes (min-max)
#SBATCH --ntasks-per-node=1 # Number of tasks per node (max)
#SBATCH --ntasks=1 # Number of tasks (processes)
#SBATCH --cpus-per-task=4 # Number of cores per task (threads)

#Descriptio: This script will run de novo genome assembly using SPAdes --rnaviral
#Usage: srun spades_viralrna.sh

echo "Running BEAST Coalescent Analysis"

module load java/17.0.3
module load beast2/2.7.3 
module load beagle/4.0.1
# Load software from Kamiak repository
srun beast -beagle_SSE -threads 12 -overwrite bdsky_D2.xml    
# Each task runs this program (total 1 times)
# Each srun is a job step, and spawns ntasks
echo "Completed job $SLURM_JOBID on nodes $SLURM_JOB_NODELIST"

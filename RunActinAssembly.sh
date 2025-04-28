#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=1GB
#SBATCH --job-name=Actin
#SBATCH --output=slurm_%j.out
#SBATCH --array=1-10
#SBATCH --account=pi-emunro

module load python
source activate /project2/emunro/MyEnv 
module unload python

cd Actin/ActinStructure/Python-Cpp
python AllActinMixedNucleates.py $SLURM_ARRAY_TASK_ID

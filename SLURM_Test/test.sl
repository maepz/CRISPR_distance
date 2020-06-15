#!/bin/bash 
#SBATCH --job-name=Testing123
#SBATCH --account=def-bacc
#SBATCH --time=0:0:30
#SBATCH --mem-per-cpu=32M
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END
#SBATCH --error='/home/maeperez/projects/def-bacc/maeperez/CRISPR_distance/SLURM/JOBNAME-'%A_%a'.err' 
#SBATCH --output='/home/maeperez/projects/def-bacc/maeperez/CRISPR_distance/SLURM/JOBNAME-'%A_%a'.out' 
#SBATCH --mail-user=maepz@hotmail.com

source /home/maeperez/virtualenv/py36/bin/activate
python test.py

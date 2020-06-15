#!/bin/bash 
#SBATCH --job-name=CRISPR_Distances
#SBATCH --account=def-bacc
#SBATCH --time=6:0:0
#SBATCH --mem-per-cpu=2048M
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END
#SBATCH --error='/home/maeperez/projects/def-bacc/maeperez/CRISPR_distance/SLURM/JOBNAME-'%A_%a'.err' 
#SBATCH --output='/home/maeperez/projects/def-bacc/maeperez/CRISPR_distance/SLURM/JOBNAME-'%A_%a'.out' 
#SBATCH --mail-user=maepz@hotmail.com

source /home/maeperez/virtualenv/py36/bin/activate
python Run_Parallel_CRISPR_Distance.py Mydata/arrays_in_whole_dataset_v6.txt 

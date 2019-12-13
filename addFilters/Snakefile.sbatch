#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH -e /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_addFilter.err
#SBATCH -o /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_addFilter.out
#SBATCH --mem=8G
#SBATCH --partition=hbfraser,owners,hns,normal
#SBATCH --time=24:00:00
#SBATCH -J s_addFilters

module load fraserconda

cd ~/scripts/FraserLab/somaticMutationsProject/addFilters/

snakemake --printshellcmds --keep-going --max-jobs-per-second 3 --max-status-checks-per-second 0.016 --nolock --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --cluster "../submit.py"


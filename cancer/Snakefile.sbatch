#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH -e /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_cancer.err
#SBATCH -o /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_cancer.out
#SBATCH --mem=4G
#SBATCH --partition=hbfraser,hns,normal,owners
#SBATCH --time=24:00:00
#SBATCH -J s_cancer

module load anaconda3
source activate fraserconda

cd ~/scripts/FraserLab/somaticMutationsProject/cancer/

snakemake --nolock --printshellcmds --keep-going --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --cluster "../submit.py"

#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH -e /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_chromatin.err
#SBATCH -o /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_chromatin.out
#SBATCH --mem=4G
#SBATCH --partition=hbfraser,hns,normal
#SBATCH --time=24:00:00
#SBATCH -J s_chromatin


module load fraserconda
source activate fraserconda

cd ~/scripts/FraserLab/somaticMutationsProject/chromatin/

snakemake --restart-times 1 --nolock --printshellcmds --keep-going --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --cluster "../submit.py"

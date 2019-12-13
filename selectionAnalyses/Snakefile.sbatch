#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH -e /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_selection.err
#SBATCH -o /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_selection.out
#SBATCH --mem=4G
#SBATCH --partition=hbfraser,owners,hns,normal
#SBATCH --time=24:00:00
#SBATCH -J s_selection


module load fraserconda
source activate fraserconda

cd ~/scripts/FraserLab/somaticMutationsProject/selectionAnalyses/

snakemake --max-jobs-per-second 3 --max-status-checks-per-second 0.016 --cluster-config ../cluster.json --cluster-status jobState --jobs 1000 --keep-going --cluster "../submit.py"

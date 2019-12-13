#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH -e /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_clusterMutAnalyses.err
#SBATCH -o /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_clusterMutAnalyses.out
#SBATCH --mem=4G
#SBATCH --partition=hbfraser
#SBATCH --time=8:00:00

module load anaconda3

cd /home/users/paedugar/scripts/FraserLab/somaticMutationsProject/tempAnalyses/clusterMutations

snakemake --nolock --keep-going --cluster-config ../../cluster.json --cluster-status jobState --jobs 500 --cluster "../../submit.py" -s Snakefile_separateClusterVsNonCluster.smk
#snakemake --nolock --keep-going --cluster-config ../../cluster.json --cluster-status jobState --jobs 1900 --cluster "../../submit.py"

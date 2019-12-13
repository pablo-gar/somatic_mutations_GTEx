#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH -e /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_seq_depth.err
#SBATCH -o /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_seq_depth.out
#SBATCH --mem=8G
#SBATCH --partition=hbfraser
#SBATCH --time=9:00:00
#SBATCH -J s_seq_depth

module load fraserconda
source activate fraserconda

cd ~/scripts/FraserLab/somaticMutationsProject/mappingMutationCalling_pipeline/

snakemake --keep-going --snakefile SnakefileSeqdepth.smk --max-jobs-per-second 3 --max-status-checks-per-second 0.016 --nolock --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --cluster "../submit.py"


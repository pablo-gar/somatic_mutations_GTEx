#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH -e /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_vaf.err
#SBATCH -o /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_vaf.out
#SBATCH --mem=4G
#SBATCH --partition=hbfraser,hns,normal
#SBATCH --time=12:00:00


####
# MODS

module load anaconda3
source activate fraserconda

cd ~/scripts/FraserLab/somaticMutationsProject/vaf/

snakemake --nolock --printshellcmds --keep-going --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --cluster "../submit.py"
#snakemake --nolock --printshellcmds --keep-going --cluster-config ../cluster.json --cluster-status jobState --jobs 7 --cluster "../submit.py"


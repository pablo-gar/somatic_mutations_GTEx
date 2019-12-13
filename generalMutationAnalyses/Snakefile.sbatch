#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH -e /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_generalAnalyses.err
#SBATCH -o /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_generalAnalyses.out
#SBATCH --mem=4G
#SBATCH --partition=hbfraser,hns,normal
#SBATCH --time=4:00:00


####
# MODS

#PATH=$PATH:$HOME/bin:$HOME/vcftools_0.1.13/bin:$HOME/bcftools_1.6/bin
PATH=$HOME/bin:$PATH:$HOME/.local/bin:$HOME/gatk-4.0.3.0:$HOME/samtools_1.6/bin
export PATH

MODULEPATH=$MODULEPATH:/share/PI/hbfraser/modules/modules
export MODULEPATH
#source /home/groups/hbfraser/shared_environment/fraser_lab_env.sh

# END MODS
######

module load anaconda3
source activate fraserconda

cd ~/scripts/FraserLab/somaticMutationsProject/generalMutationAnalyses/

#snakemake --nolock --restart-times 1 --printshellcmds --keep-going --cluster-config ../cluster.json --cluster-status jobState --jobs 7 --cluster "../submit.py"
snakemake --nolock --printshellcmds --keep-going --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --cluster "../submit.py"


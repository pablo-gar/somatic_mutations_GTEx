#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH -e /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_mutationCalling.err
#SBATCH -o /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_mutationCalling.out
#SBATCH --mem=4G
#SBATCH --partition=hbfraser,hns,normal
#SBATCH --time=72:00:00

#module load fraserconda/5.0

#####
## MODULES
PATH=$HOME/bin:$PATH:$HOME/.local/bin:$HOME/gatk-4.0.3.0:$HOME/samtools_1.6/bin
export PATH

MODULEPATH=$MODULEPATH:/share/PI/hbfraser/modules/modules
export MODULEPATH
#source /home/groups/hbfraser/shared_environment/fraser_lab_env.sh

module load fraserconda
source activate fraserconda
## END MODULES
#####

cd ~/scripts/FraserLab/somaticMutationsProject/mutationCalling/
date
echo "Start snakemake"

snakemake --keep-going --max-jobs-per-second 15 --restart-times 2 --max-status-checks-per-second 0.016 --nolock --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --cluster "../submit.py"
#while :
#do
#    snakemake --keep-going --max-jobs-per-second 3 --max-status-checks-per-second 0.016 --nolock --restart-times 2 --cluster-config ../cluster.json --cluster-status jobState --jobs 1900 --cluster "../submit.py"
#    if [ $? == 0 ]
#    then
#        break
#    fi
#    
#    date
#    echo "Something went wrong, resubmmiting snakemake"
#done

date
echo "Snakemake done!"

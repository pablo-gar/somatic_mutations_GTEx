#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH -e /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_mapping.err
#SBATCH -o /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_mapping.out
#SBATCH --mem=8G
#SBATCH --partition=hbfraser
#SBATCH --time=9:00:00
#SBATCH -J s_mapping

#module load fraserconda/5.0
module load fraserconda
source activate fraserconda

cd ~/scripts/FraserLab/somaticMutationsProject/mappingMutationCalling_pipeline/

#-------------------------
# Fastq prep
date
echo -e "Starting fastq prep\n"
# Fastq prep
snakemake --snakefile SnakefileFastqPrep.smk --max-jobs-per-second 3 --max-status-checks-per-second 0.016 --nolock --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --cluster "../submit.py"

if [ $? != 0 ]
then
    echo -e "ERROR something went wrong with fastq prep\n"
    exit 1
fi
date
echo -e "Fastq prep Snakemake done!\n\n"

##-------------------------
## Create pileups
date
echo -e "Starting Pileups snakemake\n"
#while :
#do
    snakemake --snakefile SnakefilePileup.smk --printshellcmds --keep-going --max-jobs-per-second 3 --max-status-checks-per-second 0.016 --nolock --restart-times 2 --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --cluster "../submit.py"
    #snakemake --snakefile SnakefilePileup.smk --printshellcmds --max-jobs-per-second 3 --max-status-checks-per-second 0.016 --nolock --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --cluster "../submit.py"
    
#    if [ $? == 0 ]
#    then
#        break
#    fi
#    
#    date
#    echo -e "ERROR something went wrong with pileups creation. Re-trying\n"
    
#done

date
echo -e "Create pileups Snakemake done!\n\n"


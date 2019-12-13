#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH -e /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_signatures_bayes.err
#SBATCH -o /scratch/users/paedugar/somaticMutationsProject/clusterFiles/1_Snakemake_signatures_bayes.out
#SBATCH --mem=8G
#SBATCH --partition=hbfraser
#SBATCH --time=72:00:00

module load fraserconda
source activate fraserconda

cd ~/scripts/FraserLab/somaticMutationsProject/mutationSignatures_bayes/

#snakemake  --rerun-incomplete --cluster-config cluster.json --cluster-status jobState --jobs 1900 --cluster "sbatch --parsable --partition={cluster.partition} --time={cluster.time} --cpus-per-task={cluster.cpus} --mem={cluster.mem} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error}"
#snakemake --cluster-config cluster.json --cluster-status jobState --jobs 1900 --cluster "sbatch --parsable --partition={cluster.partition} --time={cluster.time} --cpus-per-task={cluster.cpus} --mem={cluster.mem} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error}"

snakemake --keep-going --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --max-jobs-per-second 3 --max-status-checks-per-second 0.016 --nolock --cluster "../submit.py"

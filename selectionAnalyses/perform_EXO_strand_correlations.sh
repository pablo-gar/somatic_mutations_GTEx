module load R/3.4.0
cd R/

# Get dnds and annotated mutations
Rscript dndsInMutationMaps.R 4 1 0.8 /scratch/users/paedugar/somaticMutationsProject/selection/dndsout/Whole_Blood_EXO/n6_0.0_0.7/dndsout_ /scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Whole_Blood_EXO/n6_0.0_0.7/*

# Get strand fc stats
Rscript plot_strandDifferences.R /scratch/users/paedugar/somaticMutationsProject/selection/dndsout/Whole_Blood_EXO/n6_0.0_0.7/dndsout_\,geneExpression\:AllExpressed\,selectionMaf\:0%_100%\,_genesAllMutations.txt /scratch/users/paedugar/somaticMutationsProject/selection/strand_differences/Whole_Blood_EXO-n6_0.0_0.7.normalized.pdf /scratch/users/paedugar/somaticMutationsProject/selection/strand_differences/Whole_Blood_EXO-n6_0.0_0.7.raw.pdf /scratch/users/paedugar/somaticMutationsProject/selection/strand_differences/Whole_Blood_EXO-n6_0.0_0.7.box.normalized.pdf /scratch/users/paedugar/somaticMutationsProject/selection/strand_differences/Whole_Blood_EXO-n6_0.0_0.7.box.raw.pdf /scratch/users/paedugar/somaticMutationsProject/selection/strand_differences/Whole_Blood_EXO-n6_0.0_0.7.box2.raw.pdf /scratch/users/paedugar/somaticMutationsProject/selection/strand_differences/Whole_Blood_EXO-n6_0.0_0.7.raw_strand_mut_avg.txt /scratch/users/paedugar/somaticMutationsProject/selection/strand_differences/Whole_Blood_EXO-n6_0.0_0.7.age_cor_fc_strand.txt /scratch/users/paedugar/somaticMutationsProject/selection/strand_differences/Whole_Blood_EXO-n6_0.0_0.7.per_gene_fc_strand.txt

# Get correlations
mkdir -p /scratch/users/paedugar/somaticMutationsProject/selection/validation_strand_fc_perGene_blood/
Rscript compare_2_tissues_gene_fc.R /scratch/users/paedugar/somaticMutationsProject/selection/strand_differences/Whole_Blood-n6_0.0_0.7.per_gene_fc_strand.txt /scratch/users/paedugar/somaticMutationsProject/selection/strand_differences/Whole_Blood_EXO-n6_0.0_0.7.per_gene_fc_strand.txt RNA_seq exome_seq /scratch/users/paedugar/somaticMutationsProject/selection/validation_strand_fc_perGene_blood/

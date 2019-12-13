library("ggplot2")
library("dplyr")
library("tidyr")
source("../../../R/ggthemes.R")

main <- function(cmdArgs = commandArgs(T)){
    
    overlap_file <- cmdArgs[1]
    out_prefix <- cmdArgs[2]
    
    #overlap_file <- "/scratch/users/paedugar/somaticMutationsProject/comparison_yizhak/tables/table_mutations_all_filters_yizhakOverlap.bed"
    #out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/comparison_yizhak/filter_counts/"
    
    
    # Read file
    overlap <- read.delim(overlap_file, sep="\t", stringsAsFactors=F, header=F)
    overlap$tissue_mine <- TISSUES[overlap$V28]
    
    # Get only in same tissue overlaps
    overlap_same <- overlap[overlap$V17 == overlap$tissue_mine & !is.na(overlap$tissue_mine),]
    
    # Mutations that did not pass my filters
    overlap_miss <- overlap_same[!overlap_same$V27 %in% c("PASS", "clustered_mutation") ,]
    
    
    # Collect individual mutations from yizhak and see their associated filter
    
    id <- with(overlap_miss, paste0(V1, '.', V2, '.', V6, '.', V7, '.', V8))
    overlap_miss$id <- id
    
    all_filters <- strsplit(overlap_miss$V27, ';')
    uniq_muts_filters <- matrix(0, nrow = length(unique(id)), ncol = length(unique(unlist(all_filters))), dimnames = list(unique(id), unique(unlist(all_filters))))
    
    for(i in 1:length(all_filters)) {
        current_id <- overlap_miss$id[i]
        current_filters <- all_filters[[i]]
        
        uniq_muts_filters[current_id, current_filters] <- uniq_muts_filters[current_id, current_filters] + 1
    }
    
    uniq_muts_filters <- as.data.frame(uniq_muts_filters)
    uniq_muts_filters$id <- rownames(uniq_muts_filters) 
    uniq_muts_filters <- uniq_muts_filters[, colnames(uniq_muts_filters) != "clustered_mutation"]
    
    # Do plots
    filter_count_apperance <- uniq_muts_filters %>%
        gather('filter_name', 'value', -id) %>%
        group_by(filter_name) %>%
        summarise(total = sum(value > 0)) %>%
        ungroup()
    
    p <- ggplot(filter_count_apperance, aes(x=filter_name, y=total)) +
            geom_bar(stat = 'identity') +
            coord_flip() +
            theme_grid_y()
        
    ggsave(paste0(out_prefix, "number_yizhakMuts_with_filter.pdf"), p)
    write.table(uniq_muts_filters, paste0(out_prefix, "filter_counts_per_mut.txt"), sep = "\t", quote=F, row.names=F, col.names=T)
    
}

# Conversion of tissues
TISSUES = c(Adipose_Subcutaneous = "AdiposeTissue",
                      Adipose_Visceral_Omentum = "AdiposeTissue",
                      Adrenal_Gland = "AdrenalGland",
                      Artery_Aorta = "Artery_Aorta",
                      Artery_Coronary = "Artery_Coronary",
                      Artery_Tibial = "Artery_Tibial",
                      Brain_Caudate_basal_ganglia = "Brain",
                      Brain_Cortex = "Brain",
                      Brain_Frontal_Cortex_BA9 = "Brain",
                      Brain_Hippocampus = "Brain",
                      Brain_Hypothalamus = "Brain",
                      Brain_Nucleus_accumbens_basal_ganglia = "Brain",
                      Brain_Putamen_basal_ganglia = "Brain",
                      Breast_Mammary_Tissue = "Breast",
                      Colon_Sigmoid = "Colon",
                      Colon_Transverse = "Colon",
                      Esophagus_Gastroesophageal_Junction = "Esophagus",
                      Esophagus_Mucosa = "Esophagus",
                      Esophagus_Muscularis = "Esophagus",
                      Heart_Atrial_Appendage = "Heart",
                      Heart_Left_Ventricle = "Heart",
                      Liver = "Liver",
                      Lung = "Lung",
                      Muscle_Skeletal = "Muscle",
                      Nerve_Tibial = "Nerve",
                      Ovary = "Ovary",
                      Pancreas = "Pancreas",
                      Pituitary = "Pituitary",
                      Prostate = "Prostate",
                      Skin_Not_Sun_Exposed_Suprapubic = "Skin",
                      Skin_Sun_Exposed_Lower_leg = "Skin",
                      Small_Intestine_Terminal_Ileum = "SmallIntestine",
                      Spleen = "Spleen",
                      Stomach = "Stomach",
                      Thyroid = "Thyroid",
                      Whole_Blood = "BloodVessel")

main()

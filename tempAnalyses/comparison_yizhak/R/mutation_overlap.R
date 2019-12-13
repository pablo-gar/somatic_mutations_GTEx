# Calculates level of overlap with yizhak mutations
library("dplyr")
library("VennDiagram")
library("ggplot2")
source("../../../R/ggthemes.R")

main <- function(cmdArgs = commandArgs(T)) {
    
    garcia_mut_file <- cmdArgs[1]
    yizhak_mut_file <- cmdArgs[2]
    
    garcia_alt_col <- as.numeric(cmdArgs[3])
    yizhak_alt_col <- as.numeric(cmdArgs[4])
    garcia_tis_col <- as.numeric(cmdArgs[5])
    yizhak_tis_col <- as.numeric(cmdArgs[6])
    
    out_prefix <- cmdArgs[7]
    

    #DEBUG ARGS
    #garcia_mut_file <- "/scratch/users/paedugar/somaticMutationsProject/comparison_yizhak/tables/table_mutations_all.bed"
    #yizhak_mut_file <- "/scratch/users/paedugar/somaticMutationsProject/comparison_yizhak/tables/table_mutations_yizhak.bed"
    #garcia_alt_col <- 5
    #yizhak_alt_col <- 7
    #garcia_tis_col <- 9
    #yizhak_tis_col <- 17
    

    
    
    # Read tables
    garcia_mut <- read.delim(garcia_mut_file, sep = "\t", header = F, stringsAsFactors = F)
    yizhak_mut <- read.delim(yizhak_mut_file, sep = "\t", header = F, stringsAsFactors = F)
    
    garcia_ncol <- ncol(garcia_mut)
    yizhak_ncol <- ncol(yizhak_mut)
    
    # Get overlap
    mut_overlap <- get_overlap(garcia_mut, yizhak_mut)
    
    # Creates venn diagrams
    pdf(paste0(out_prefix, "overlap_all_mutations.pdf"), width = 4, height = 4)
    make_venn(mut_overlap, garcia_ncol, yizhak_ncol, garcia_alt_col, yizhak_alt_col, unique_muts = F)
    dev.off()
    
    pdf(paste0(out_prefix, "overlap_unique_mutations.pdf"), width = 4, height = 4)
    make_venn(mut_overlap, garcia_ncol, yizhak_ncol, garcia_alt_col, yizhak_alt_col, unique_muts = T)
    dev.off()
    
    # Gets overlap percentage per tissue
    perTissue_overlap <- get_perTissue_overlap(mut_overlap, 
                                               garcia_ncol, yizhak_ncol, 
                                               garcia_alt_col, yizhak_alt_col, 
                                               garcia_tis_col, yizhak_tis_col, 
                                               unique_muts = F)
    
    perTissue_overlap_unique <- get_perTissue_overlap(mut_overlap, 
                                               garcia_ncol, yizhak_ncol, 
                                               garcia_alt_col, yizhak_alt_col, 
                                               garcia_tis_col, yizhak_tis_col, 
                                               unique_muts = T)
    
    # Plot percentage
    p1 <- ggplot(perTissue_overlap, aes(x = tissue, y = percent)) +
        geom_bar(stat = "identity", colour = "black", fill = "lightblue") +
        geom_text(aes(label = label), hjust = 0.5, vjust = 0) +
        labs(subtitle = "n is the total number of mutations in Yizhak et al.", x = "", y = "Percent of Yizhak et al. mutations in Garcia-Nieto et al.") +
        theme_grid_y() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1)) 
    
    p2 <- ggplot(perTissue_overlap_unique, aes(x = tissue, y = percent)) +
        geom_bar(stat = "identity", colour = "black", fill = "lightblue") +
        geom_text(aes(label = label), hjust = 0.5, vjust = 0) +
        labs(subtitle = "n is the total number of mutations in Yizhak et al.", x = "", y = "Percent of Yizhak et al. mutations in Garcia-Nieto et al.") +
        theme_grid_y() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1)) 
    
    ggsave(paste0(out_prefix, "perTissue_percent_overlap_all.pdf"), p1, width = 12, height = 6)
    ggsave(paste0(out_prefix, "perTissue_percent_overlap_unique.pdf"), p2, width = 12, height = 6)
    
}

get_overlap <- function(x, y) {
    
    x$id <- paste(x[,1], x[,2])
    y$id <- paste(y[,1], y[,2])
    
    result <- full_join(x, y, by = "id")
    result <- select(result, -id)
    
    return(result)
    
}

get_overlap_counts <- function(xy, x.ncol, y.ncol, x.alt_col, y.alt_col, unique_muts = F) {
    x.col <- 1
    y.col <- x.ncol + 1
    y.alt_col <- x.ncol + y.alt_col
    
    if(unique_muts) {
        xy <- xy[!duplicated(paste(xy[,x.col], xy[,x.col + 1], xy[,x.alt_col])) | 
                 is.na(xy[,x.col]), ]
             
        xy <- xy[!duplicated(paste(xy[,y.col], xy[,y.col + 1], xy[,y.alt_col])) |
                 is.na(xy[,y.col]), ]
    }
    
    x.total <- !is.na(xy[,x.col])
    y.total <- !is.na(xy[,y.col])
    overlap_count <- x.total & y.total
    
    return(c(x = sum(x.total), y = sum(y.total), overlap = sum(overlap_count)))
    
}

make_venn <- function(xy, x.ncol, y.ncol, x.alt_col, y.alt_col, unique_muts = F) {
    
    overlap_counts <- get_overlap_counts(xy, x.ncol, y.ncol, x.alt_col, y.alt_col, unique_muts)
    
    draw.pairwise.venn(overlap_counts["x"], overlap_counts["y"], overlap_counts["overlap"], category = c("Garcia-Nieto et al.", "Yizhak et al."), 
                       fill = c("lightsteelblue3", "tomato3"), lwd = rep(0.5,2),
                       alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
    
}

get_perTissue_overlap <- function(xy, x.ncol, y.ncol, x.alt_col, y.alt_col, x.tis_col, y.tis_col,  unique_muts = F) {
    
    # Conversion of tissues
    tissue_conversion = c(Adipose_Subcutaneous = "AdiposeTissue",
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
    
    x.col <- 1
    y.col <- x.ncol + 1
    y.alt_col <- x.ncol + y.alt_col
    y.tis_col <- x.ncol + y.tis_col
    
    if(unique_muts) {
        xy <- xy[!duplicated(paste(xy[,x.col], xy[,x.col + 1], xy[,x.alt_col])) | 
                 is.na(xy[,x.col]), ]
             
        xy <- xy[!duplicated(paste(xy[,y.col], xy[,y.col + 1], xy[,y.alt_col])) |
                 is.na(xy[,y.col]), ]
    }
    
    # Getting shared tissues
    xy[,x.tis_col] <- tissue_conversion[xy[,x.tis_col]]
    shared_tis <-  unique(xy[ ,x.tis_col]) [unique(xy[ ,x.tis_col]) %in% unique(xy[ ,y.tis_col])]
    shared_tis <- shared_tis[!is.na(shared_tis)]
    results <- list()
    for(tis in unique(shared_tis)) {
        
        both <- xy[,x.tis_col] ==  xy[,y.tis_col] & xy[,x.tis_col] == tis
        only_x <- xy[,x.tis_col] == tis &  is.na(xy[,y.tis_col])
        only_y <- xy[,y.tis_col] == tis &  is.na(xy[,x.tis_col])
        
        both[is.na(both)] <- F
        
        current_xy <- xy[ both | only_x | only_y, ]
    
        x.total <- !is.na(current_xy[,x.col])
        y.total <- !is.na(current_xy[,y.col])
        overlap_count <- x.total & y.total
        
        results[[tis]] <- data.frame(x = sum(x.total), y = sum(y.total), overlap = sum(overlap_count), tissue = tis, stringsAsFactors = F)
    }
    
    
    results <- do.call(rbind, results)
    
    results$percent <- results$overlap / results$y * 100 
    results$label <- paste("n =", results$y)
    results$tissue <- factor(results$tissue, levels = results$tissue[order(-results$percent)], ordered = T)
    
    return(results)
    
}

main()

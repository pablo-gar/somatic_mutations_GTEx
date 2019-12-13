# Performs basic analysis of VAF and metadata
source("../../R/plots.R", chdir=T)
source("../../R/ggthemes.R", chdir=T)
library("dplyr")
library("tidyr")
library("ggplot2")

main <- function(cmdArgs = commandArgs(T)) {
    
    mutation_map_file <- cmdArgs[1]
    out_folder <- cmdArgs[2]
    
    mutation_map_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/1_final_tables/all_mutations_annotated.txt"
    metadata_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/1_final_tables/metadata.txt"
    out_folder <- "/scratch/users/paedugar/somaticMutationsProject/vaf/general_vaf_analysis/"
    
    ## READ mutations
    mutation_map <- read.table(mutation_map_file, sep="\t", stringsAsFactors=F, header=T)
    
    ## READ metadata
    metadata <- read.table(file=metadata_file, header=T, stringsAsFactors=F)
    
    
    
    #---------------------------------------# 
    # Simple VAF analysis   
    # Get multiple organizations of the data
    
    # median, mean per sample
    vaf_per_sample <- get_vaf_per_sample(mutation_map, vaf, metadata)
    
    # median, mean per sample vaf corrected 
    vaf_corrected_per_sample <- get_vaf_per_sample(mutation_map, vaf_corrected, metadata)
    
    
    # median, mean per sample per mutation type 
    vaf_per_sample_per_mut <- get_vaf_per_sample_per_mut(mutation_map, vaf, metadata)    
    vaf_corrected_per_sample_per_mut <- get_vaf_per_sample_per_mut(mutation_map, vaf_corrected, metadata)
    
    # median, mean per tissue
    vaf_per_tissue <- mutation_map %>%
        group_by(tissue) %>%
        summarise(vaf_median = median(vaf), vaf_mean = mean(vaf)) %>%
        ungroup()
    #---------------------------------------# 
    
    
    
    # PLOTS 
    
    #---------------------------------------# 
    # Plot VAFs vs coverage
    p_vaf_coverage <- scatter(mutation_map, x="coverage", y="vaf", alpha=0.05) + theme_sleek()
    p_vaf_corrected_coverage <- scatter(mutation_map, x="coverage", y="vaf_corrected", alpha=0.05) + theme_sleek()
    p_vaf_coverage_per_mut <- scatter(mutation_map, x="coverage", y="vaf", facet_x="mut", alpha=0.05) + theme_sleek()
    
    p_vaf_coverage_density <- ggplot(mutation_map, aes(x=coverage, y=vaf)) +
        stat_binhex(bins=100, aes(fill=log10(..count..))) +
        theme_sleek()
    
    p_vaf_coverage_density <- ggplot(mutation_map, aes(x=coverage, y=vaf)) +
        stat_binhex(bins=100) +
        theme_sleek()
    
    p_vaf_coverage_per_mut_density <- ggplot(mutation_map, aes(x=coverage, y=vaf)) +
        stat_binhex() +
        facet_grid(.~mut) +
        theme_sleek()
    
    ggsave(file.path(out_folder, "vaf_vs_coverage.png"), p_vaf_coverage )
    ggsave(file.path(out_folder, "vaf_corrected_vs_coverage.png"), p_vaf_corrected_coverage)
    ggsave(file.path(out_folder, "vaf_vs_coverage_per_mut.png"), p_vaf_coverage_per_mut, width = 14)
    ggsave(file.path(out_folder, "vaf_vs_coverage_density.pdf"), p_vaf_coverage_density) 
    #---------------------------------------# 
    
    
        
    #---------------------------------------# 
    # Plot median VAF per sample and transcriptome diversity
    p_per_sample_vaf_seqDepth <- plot_vaf_vs_seqDepth_transDiversity(vaf_per_sample)
    p_per_sample_vaf_corrected_seqDepth <- plot_vaf_vs_seqDepth_transDiversity(vaf_corrected_per_sample)
    
    ggsave(file.path(out_folder, "vaf_per_sample_vs_seq_depth.pdf"), p_per_sample_vaf_seqDepth, height=9, width=9)
    ggsave(file.path(out_folder, "vaf_corrected_per_sample_vs_seq_depth.pdf"), p_per_sample_vaf_corrected_seqDepth, height=9, width=9)
    #---------------------------------------# 
                                     
                                     
    #---------------------------------------# 
    # Plot median VAF per tissue
    p_vaf_median_per_tissue <- plot_vaf_dist_per_tissue(vaf_per_sample, vaf_median)
    p_vaf_mean_per_tissue <- plot_vaf_dist_per_tissue(vaf_per_sample, vaf_mean)
    p_vaf_corrected_median_per_tissue <- plot_vaf_dist_per_tissue(vaf_corrected_per_sample, vaf_median)
    p_vaf_corrected_mean_per_tissue <- plot_vaf_dist_per_tissue(vaf_corrected_per_sample, vaf_mean)
    
    ggsave(file.path(out_folder, "vaf_per_tissue_median.pdf"), p_vaf_median_per_tissue, height=7, width=12)
    ggsave(file.path(out_folder, "vaf_per_tissue_mean.pdf"), p_vaf_mean_per_tissue, height=7, width=12)
    ggsave(file.path(out_folder, "vaf_corrected_per_tissue_median.pdf"), p_vaf_corrected_mean_per_tissue, height=7, width=12)
    ggsave(file.path(out_folder, "vaf_corrected_per_tissue_mean.pdf"), p_vaf_corrected_mean_per_tissue, height=7, width=12)
    #---------------------------------------# 
        
    
    #---------------------------------------# 
    # Plot features associations per tissue
    
    p_sample_metadata_per_tissue_cors <- plot_cors_metadata(vaf_per_sample)
    p_sample_metadata_per_tissue_cors_corrected <- plot_cors_metadata(vaf_corrected_per_sample)
    
    ggsave(file.path(out_folder, "vaf_cor_metadata.pdf"), p_sample_metadata_per_tissue_cors, width = 14, height=7 )
    ggsave(file.path(out_folder, "vaf_corrected_cor_metadata.pdf"), p_sample_metadata_per_tissue_cors_corrected, width = 14, height=7 )
    
    #---------------------------------------# 
    
    
    #---------------------------------------# 
    # Plot features associations per tissue per mutation
    p_sample_metadata_per_tissue_cors_per_mut <- plot_cors_metadata_per_mut(vaf_per_sample_per_mut)
    p_sample_metadata_per_tissue_cors_per_mut_corrected <- plot_cors_metadata_per_mut(vaf_corrected_per_sample_per_mut)
    
    for(i in names(p_sample_metadata_per_tissue_cors_per_mut))
        ggsave(file.path(out_folder, paste0("vaf_cor_metadata_", i, ".pdf")), p_sample_metadata_per_tissue_cors_per_mut[[i]], width = 14, height=7)
    
    for(i in names(p_sample_metadata_per_tissue_cors_per_mut))
        ggsave(file.path(out_folder, paste0("vaf_corrected_cor_metadata_", i, ".pdf")), p_sample_metadata_per_tissue_cors_per_mut_corrected[[i]], width = 14, height=7)
    #---------------------------------------# 
    
    
}

# Plot features associations per tissue
plot_cors_metadata <- function(x) {
        
    toPlot <- gather(x, "feature", "value", n_uniqueMapped, transcriptome_diversity, AGE, GENDER, BMI, C1, C2)
    spearman_cors <- toPlot %>%
    group_by(tissue, feature) %>%
    summarise(spearman=cor.test(value, vaf_median, method="spearman")[["estimate"]], 
              pvalue=cor.test(value, vaf_median, method="spearman")[["p.value"]],
              log10_signed_pvalue = -log10(pvalue) *(abs(spearman)/ spearman)) %>%
    ungroup()
    
    ggplot(spearman_cors, aes(x=log10_signed_pvalue, y=tissue)) +
        geom_point() +
        geom_vline(xintercept=0, colour="grey50", linetype='dashed') + 
        facet_grid(.~feature) +
        theme_sleek() +
        theme_grid()
        
}


# Get vaf per sample
get_vaf_per_sample <- function(x, metric, metadata) {
    
    metric <- enquo(metric)
    
    vaf_per_sample <- x %>%
        group_by(sraIds) %>%
        summarise(gtexIds = gtexIds[1], gtexIds_samples = gtexIds_samples[1], tissue = tissue[1],
                  vaf_median = median(!!metric), vaf_mean = mean(!!metric)) %>%
        ungroup()
    
    vaf_per_sample <- left_join(vaf_per_sample, metadata, by="sraIds")
    
    return(vaf_per_sample)
    
}

# Get vaf per sample per mutation type 
get_vaf_per_sample_per_mut <- function(x, metric, metadata) {
    
    metric <- enquo(metric)
    
    vaf_per_sample_per_mut <- x %>%
        group_by(sraIds, mut) %>%
        summarise(gtexIds = gtexIds[1], gtexIds_samples = gtexIds_samples[1], tissue = tissue[1],
                        vaf_median = median(!!metric), vaf_mean = mean(!!metric)) %>%
        ungroup()

    vaf_per_sample_per_mut <- left_join(vaf_per_sample_per_mut, metadata, by="sraIds")
    
    return(vaf_per_sample_per_mut)
}


# Plot features associations per tissue per mutation type
plot_cors_metadata_per_mut <- function(x) {
    
    toPlot <- gather(x, "feature", "value", n_uniqueMapped, transcriptome_diversity, AGE, GENDER, BMI, C1, C2)
    spearman_cors <- toPlot %>%
        group_by(tissue, feature, mut) %>%
        summarise(spearman=cor.test(value, vaf_median, method="spearman")[["estimate"]], 
                  pvalue=cor.test(value, vaf_median, method="spearman")[["p.value"]],
                  log10_signed_pvalue = -log10(pvalue) *(abs(spearman)/ spearman)) %>%
        ungroup()
    
    p_sample_metadata_per_tissue_cors_per_mut <- list()
    for(i in unique(toPlot$mut)) {
        current <- spearman_cors[spearman_cors$mut == i, ]
        p_sample_metadata_per_tissue_cors_per_mut[[i]] <- ggplot(current, aes(x=log10_signed_pvalue, y=tissue)) +
            geom_point() +
            facet_grid(.~feature) +
            ggtitle(i) +
            theme_sleek() +
            theme_grid()
        
        #dev.new()
        #print(p_sample_metadata_per_tissue_cors_per_mut[[i]])
    }
    
    return(p_sample_metadata_per_tissue_cors_per_mut)
}

# Plot median VAF per sample and transcriptome diversity
plot_vaf_vs_seqDepth_transDiversity <- function(x) {

    toPlot <- gather(x, "stat", "vaf", vaf_mean, vaf_median)
    toPlot <- gather(toPlot, "feature", "value", n_uniqueMapped, transcriptome_diversity)
                     
    scatter(as.data.frame(toPlot), x="value", y="vaf", scales="free", facet_x="stat", facet_y="feature", alpha=0.3) +
        theme_sleek()
}

# Plot median VAF per tissue

plot_vaf_dist_per_tissue <- function(x, metric) {
    metric <- enquo(metric)
    
    toPlot <- x
    toPlot <- toPlot %>%
        mutate(tissue = factor(tissue, ordered=T, levels=names(sort(tapply(!!metric, tissue, median)))))
                     
    ggplot(toPlot, aes(x=tissue, y=!!metric)) +
        geom_boxplot() +
        theme_sleek() +
        theme_grid_y() +
        theme(axis.text.x=element_text(angle=30, hjust=1))
}

main()

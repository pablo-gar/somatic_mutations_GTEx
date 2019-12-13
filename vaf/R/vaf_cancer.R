# Performs basic analysis of VAF and metadata
source("../../R/plots.R", chdir=T)
source("../../R/ggthemes.R", chdir=T)
source("../../R/gtex.R", chdir=T)
source("../../R/misc.R", chdir=T)
library("dplyr")
library("tidyr")
library("ggplot2")
library('ggrepel')

main <- function(cmdArgs = commandArgs(T)) {
    
    barCols <- c("#FDE725", "#440154")
    
    mutation_map_file <- cmdArgs[1]
    out_folder <- cmdArgs[2]
    
    
    mutation_map_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/1_final_tables/all_mutations_annotated.txt"
    metric <- 'vaf'
    metadata_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/1_final_tables/metadata.txt"
    out_folder <- file.path("/scratch/users/paedugar/somaticMutationsProject/", metric, "/cancer_analysis/")
    
    dir.create(out_folder, recursive=T, showWarnings=F)
    
    ## READ mutations
    mutation_map <- read.table(mutation_map_file, sep="\t", stringsAsFactors=F, header=T)
    mutation_map <- dplyr::filter(mutation_map, !is.na(cancerous) & impact %in% c('Synonymous', 'Missense', 'Nonsense'))
    
    ## READ metadata
    metadata <- read.table(file=metadata_file, header=T, stringsAsFactors=F)
    metadata$tissue <- sraToTissues(metadata$sraIds)
    
    
    ## Get metrics per tissue and sample
    vaf_cancer_per_tissue <- get_vaf_cancer_per_tissue(mutation_map, metric, metadata)
    vaf_cancer_per_tissue_gathered <- get_vaf_cancer_per_tissue_gathered(mutation_map, metric, metadata)
    vaf_cancer_per_sample <- get_vaf_cancer_per_sample(mutation_map, metric, select(metadata,-tissue))
    
    # vaf and coverage across groups of mutations based on coverage 
    vaf_cancer_coverage_groups <- get_vaf_cancer_per_group_coverage(mutation_map, metric, groups=20)
    
    # coverage and vaf fold-changes (seen in cancer / non in cancer)
    vaf_coverage_cancer_per_sample <- get_vaf_coverage_cancer_per_sample(mutation_map, metric)
    
    # Only coverage per sample
    coverage_cancer_per_sample <- get_coverage_per_sample_impact_cancer(mutation_map)
    
    vaf_cancer_per_sample_spread <- vaf_cancer_per_sample %>%
        select(-Not_in_cancer, -Seen_in_cancer) %>%
        spread(impact, fc_vaf_cancer)
    
    #---------------------------------------------
    # vaf differences acording to impact to protein and cancer status
    #---------------------------------------------
    toPlot <- mutation_map %>%
        group_by(impact, cancerous) %>%
        summarise(vaf_median = median(vaf), 
                  lower_conf=bootstrap_confidence_interval(vaf, median)[1],
                  upper_conf=bootstrap_confidence_interval(vaf, median)[2]) %>%
        ungroup()
    
    p <- ggplot(toPlot, aes(x=impact, y=vaf_median, fill=cancerous)) +
        geom_bar(colour='black', stat='identity', position=position_dodge(width=0.9)) +
        scale_fill_manual(values=barCols) +
        geom_errorbar(aes(ymax=upper_conf, ymin=lower_conf), width = 0.3, position=position_dodge(width=0.9)) +
        theme_sleek() +
        theme(axis.text.x=element_text(angle=30, hjust=1), legend.position='top')
    
    ggsave(file.path(out_folder, 'impact_cancer_vaf.pdf')) 
    
    #---------------------------------------------
    # coverage differences acording to impact to protein and cancer status
    #---------------------------------------------
    toPlot <- mutation_map %>%
        group_by(impact, cancerous) %>%
        summarise(coverage_median = median(coverage), 
                  lower_conf=bootstrap_confidence_interval(coverage, median)[1],
                  upper_conf=bootstrap_confidence_interval(coverage, median)[2]) %>%
        ungroup()
    
    p <- ggplot(toPlot, aes(x=impact, y=coverage_median, fill=cancerous)) +
        geom_bar(colour='black', stat='identity', position=position_dodge(width=0.9)) +
        scale_fill_manual(values=barCols) +
        geom_errorbar(aes(ymax=upper_conf, ymin=lower_conf), width = 0.3, position=position_dodge(width=0.9)) +
        theme_sleek() +
        theme(axis.text.x=element_text(angle=30, hjust=1), legend.position='top')
    
    ggsave(file.path(out_folder, 'impact_cancer_coverage.pdf')) 
    
    #---------------------------------------------
    # vaf differences within coverage groups
    #---------------------------------------------
    
    p <- ggplot(vaf_cancer_coverage_groups, aes(x=impact, y=vaf_median, fill=cancerous)) +
        geom_bar(stat='identity', colour='black', position=position_dodge(width=0.9)) +
        scale_fill_manual(values=barCols) +
        facet_grid(~group) +
        theme_sleek() +
        theme(axis.text.x=element_text(angle=30, hjust=1), legend.position='top')
    
    ggsave(file.path(out_folder, 'impact_cancer_vaf_binned.pdf'), width = 15) 
    
    p <- ggplot(vaf_cancer_coverage_groups, aes(x=impact, y=coverage_median, fill=cancerous)) +
        geom_bar(stat='identity', colour='black', position=position_dodge(width=0.9)) +
        scale_fill_manual(values=barCols) +
        facet_grid(~group) +
        theme_sleek() +
        theme(axis.text.x=element_text(angle=30, hjust=1), legend.position='top')
    
    ggsave(file.path(out_folder, 'impact_cancer_coverage_binned.pdf'), width=15) 
    
    #---------------------------------------------
    # Plot tissue vaf fold-chanage (cancer/non-cancer) vs seq depth
    #---------------------------------------------
    
    p <- ggplot(vaf_cancer_per_tissue, aes(x=n_uniqueMapped, y=fc_vaf_cancer)) +
        geom_point() +
        geom_text_repel(aes(label=tissue)) +
        facet_grid(.~impact) +
        geom_hline(yintercept=0, linetype='dashed') +
        theme_sleek()
    
    ggsave(file.path(out_folder, 'vaf_fc_tissues_vs_seqDepth.pdf'), width=21)
    
    #---------------------------------------------
    # cancer vs non-cancer coverage differences across tissues
    #---------------------------------------------
    
    p <- ggplot(coverage_cancer_per_sample, aes(x=impact, y=coverage_median, fill=cancerous)) +
        geom_bar(stat='identity', position=position_dodge(width=0.9)) +
        geom_errorbar(aes(ymin=lower_conf, ymax=upper_conf), width=0.4, position=position_dodge(width=0.9))+
        facet_wrap(~tissue, scale='free', ncol=6) +
        theme_sleek()
    
    ggsave(file.path(out_folder, 'cancer_non_cancer_coverage_diff_tissue.pdf'), width=21)
    
    
    #---------------------------------------------
    # cancer vs non-cancer vaf differences across tissues
    #---------------------------------------------
    
    p <- ggplot(vaf_cancer_per_tissue_gathered, aes(x=impact, y=vaf_median, fill=cancerous)) +
        geom_bar(colour='black', stat='identity', position=position_dodge(width=0.9)) +
        scale_fill_manual(values=barCols) +
        geom_errorbar(aes(ymin=lower_conf, ymax=upper_conf), width=0.4, position=position_dodge(width=0.9))+
        facet_wrap(~tissue, scale='free', ncol=6) +
        theme_sleek()
    
    ggsave(file.path(out_folder, 'cancer_non_cancer_vaf_diff_tissue.pdf'), width=21)
    
    
    #---------------------------------------# 
    # Basic correaltions
    #---------------------------------------# 
    
    
    #Syn vs Mis
    p <- scatter(as.data.frame(vaf_cancer_per_sample_spread), x='Synonymous', y='Missense', pSize=1.5, alpha=0.3) +
        theme_sleek() +
        xlab("log2(cancer synonymous VAF / non-cancer synonymous VAF)") + ylab("log2(cancer missense VAF / missense synonymous VAF)") + 
        theme_grid()
    
    ggsave(paste0(out_folder, "synonymous_vs_missense_per_sample.pdf"), p)
    
    #Non vs Mis
    p <- scatter(as.data.frame(filter(vaf_cancer_per_sample_spread, !is.na(Nonsense))), x='Synonymous', y='Missense', pSize=1.5, alpha=0.3) +
        theme_sleek() +
        xlab("log2(cancer nonsense VAF / non-cancer nonsense VAF)") + ylab("log2(cancer missense VAF / missense synonymous VAF)") + 
        theme_grid()
    
    ggsave(paste0(out_folder, "synonymous_vs_nonsense_per_sample.pdf"), p)
    
    #---------------------------------------# 
    # Metadata correlations
    #---------------------------------------# 
    
    p <- plot_cors_metadata(vaf_cancer_per_sample_spread, Missense)
    ggsave(file.path(out_folder, 'vaf_missense_cancer_metadata_cor.pdf'), width=16)
    p <- plot_cors_metadata(vaf_cancer_per_sample_spread, Synonymous)
    ggsave(file.path(out_folder, 'vaf_synonymous_cancer_metadata_cor.pdf'), width=16)
    
}

get_vaf_cancer_per_tissue_gathered <- function(x, metric, metadata) {
    
    metric <- as.name(metric)
   
    x <- x %>%
        filter(!is.na(impact)) %>%
        group_by(tissue, impact, cancerous) %>%
        summarise(vaf_median = median(!!metric),
                  lower_conf=bootstrap_confidence_interval(!!metric, median)[1],
                  upper_conf=bootstrap_confidence_interval(!!metric, median)[2],
                  ) %>%
        ungroup() 
    
    return(x)
}


get_vaf_cancer_per_group_coverage <- function(x, metric, groups = 20, by_quantile=T) {
    
    if(by_quantile) {
        groups <- quantile(x$coverage, seq(0,1, length.out=groups))
    } else {
        groups <- seq(min(x$coverage), max(x$coverage), length.out=groups)
    }
    
    x$group <- findInterval(x$coverage, groups, all.inside=T)
    
    metric <- as.name(metric)
   
    y <- x %>%
        filter(!is.na(impact)) %>%
        group_by(impact, cancerous, group) %>%
        summarise(vaf_median = median(!!metric)) %>%
        ungroup() %>%
        mutate(id = paste(impact, group, cancerous))
    
    z <- x %>%
        filter(!is.na(impact)) %>%
        group_by(impact, cancerous, group) %>%
        summarise(coverage_median = median(coverage)) %>%
        ungroup() %>%
        mutate(id = paste(impact, group, cancerous)) %>%
        select(-impact, -group, -cancerous)
            
    
    x <- left_join(y, z, by='id')
    
    
    
    
    return(x)
}

get_vaf_cancer_per_tissue <- function(x, metric, metadata) {
    
    metric <- as.name(metric)
    
    do_log <- T
    if(any(select(x,!!metric)[,1,drop=T] < 0))
        do_log <- F
    
    x <- x %>%
        filter(!is.na(impact)) %>%
        group_by(tissue, impact, cancerous) %>%
        summarise(vaf_median = median(!!metric)) %>%
        ungroup() %>%
        group_by(tissue, impact) %>%
        spread(cancerous, vaf_median) %>%
        mutate(fc_vaf_cancer = ifelse(do_log, log2(Seen_in_cancer/Not_in_cancer), Seen_in_cancer - Not_in_cancer)) %>%
        ungroup()
    
    meta <- metadata %>%
        group_by(tissue) %>%
        summarise(n_uniqueMapped = sum(as.numeric(n_uniqueMapped)), transcriptome_diversity = median(transcriptome_diversity[!is.na(transcriptome_diversity)])) %>%
        ungroup()
    
    x <- left_join(x, meta, by='tissue')
    
    
    return(x)
}

get_vaf_cancer_per_tissue <- function(x, metric, metadata) {
    
    metric <- as.name(metric)
    
    do_log <- T
    if(any(select(x,!!metric)[,1,drop=T] < 0))
        do_log <- F
   
    x <- x %>%
        filter(!is.na(impact)) %>%
        group_by(tissue, impact, cancerous) %>%
        summarise(vaf_median = median(!!metric)) %>%
        ungroup() %>%
        group_by(tissue, impact) %>%
        spread(cancerous, vaf_median) %>%
        mutate(fc_vaf_cancer = ifelse(do_log, log2(Seen_in_cancer/Not_in_cancer), Seen_in_cancer - Not_in_cancer)) %>%
        ungroup()
    
    meta <- metadata %>%
        group_by(tissue) %>%
        summarise(n_uniqueMapped = sum(as.numeric(n_uniqueMapped)), transcriptome_diversity = median(transcriptome_diversity[!is.na(transcriptome_diversity)])) %>%
        ungroup()
    
    x <- left_join(x, meta, by='tissue')
    
    
    return(x)
}

get_vaf_cancer_per_sample <- function(x, metric, metadata) {
    
    metric <- as.name(metric)
    
    do_log <- T
    if(any(select(x,!!metric)[,1,drop=T] < 0))
        do_log <- F
   
    x <- x %>%
        filter(!is.na(impact)) %>%
        group_by(impact, cancerous, sraIds) %>%
        summarise(gtexIds=gtexIds[1], tissue=tissue[1], vaf_median = median(!!metric)) %>%
        ungroup() %>%
        group_by(sraIds, impact) %>%
        spread(cancerous, vaf_median) %>%
        mutate(fc_vaf_cancer = ifelse(do_log, log2(Seen_in_cancer/Not_in_cancer), Seen_in_cancer - Not_in_cancer)) %>%
        ungroup()
    
    x <- left_join(x, metadata, by='sraIds')
    
    
    return(x)
}

get_vaf_coverage_cancer_per_sample <- function(x, metric) {
    
    metric <- as.name(metric)
    
    do_log <- T
    if(any(select(x,!!metric)[,1,drop=T] < 0))
        do_log <- F
   
    y <- x %>%
        filter(!is.na(impact)) %>%
        group_by(impact, cancerous, sraIds) %>%
        summarise(gtexIds=gtexIds[1], tissue=tissue[1], vaf_median = median(!!metric)) %>%
        ungroup() %>%
        group_by(sraIds, impact) %>%
        spread(cancerous, vaf_median) %>%
        mutate(fc_vaf_cancer = ifelse(do_log, log2(Seen_in_cancer/Not_in_cancer), Seen_in_cancer - Not_in_cancer)) %>%
        ungroup() %>%
        mutate(id = paste(impact, sraIds))
    
    z <- x %>%
        filter(!is.na(impact)) %>%
        group_by(impact, cancerous, sraIds) %>%
        summarise(coverage_median = median(coverage)) %>%
        ungroup() %>%
        group_by(sraIds, impact) %>%
        spread(cancerous, coverage_median) %>%
        mutate(fc_vaf_cancer = ifelse(do_log, log2(Seen_in_cancer/Not_in_cancer), Seen_in_cancer - Not_in_cancer)) %>%
        ungroup() %>%
        mutate(id = paste(impact, sraIds)) %>%
        select(-impact, -sraIds)
            
    
    x <- left_join(y, z, by='id')
    
    
    return(x)
}


plot_cors_metadata <- function(x, metric) {
    
    metric <- enquo(metric)
        
    toPlot <- gather(x, "feature", "value", n_uniqueMapped, transcriptome_diversity, AGE, GENDER, BMI, C1, C2)
    spearman_cors <- toPlot %>%
    group_by(tissue, feature) %>%
    summarise(spearman=cor.test(value, !!metric, method="spearman")[["estimate"]], 
              pvalue=cor.test(value, !!metric, method="spearman")[["p.value"]],
              log10_signed_pvalue = -log10(pvalue) *(abs(spearman)/ spearman)) %>%
    ungroup() %>%
    mutate(tissue=factor(tissue, ordered=T, levels=names(sort(tapply(log10_signed_pvalue[feature=="AGE"], tissue[feature=="AGE"], median)))))
    
    ggplot(spearman_cors, aes(x=log10_signed_pvalue, y=tissue)) +
        geom_point() +
        geom_vline(xintercept=0, colour="grey50", linetype='dashed') + 
        facet_grid(.~feature) +
        theme_bw() + 
        theme_sleek() + 
        theme_grid()
        
}

get_coverage_per_sample_impact_cancer <- function(x) {
    
   result <- x %>%
       group_by(impact, cancerous, tissue) %>%
       summarise(coverage_median=median(coverage),
                 lower_conf=bootstrap_confidence_interval(coverage, median)[1],
                 upper_conf=bootstrap_confidence_interval(coverage, median)[2],
                 ) %>%
       ungroup()
   
   return(result)
}

spread_n <-  function(df, key, value) {
    
    # quote key
    keyq <- rlang::enquo(key)
    # break value vector into quotes
    valueq <- rlang::enquo(value)
    s <- rlang::quos(!!valueq)
    df %>% gather(variable, value, !!!s) %>%
        unite(temp, !!keyq, variable) %>%
        spread(temp, value)
    
}

main()

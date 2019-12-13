#  Get median vaf per sample for given metric (vaf vs vaf_corrected)
get_vaf_per_sample <- function(x, metric) {
    
    x_all <- x %>%
        group_by(sraIds) %>%
        summarise(all = median(!!as.name(metric))) %>%
        ungroup() 
    
    x <- x %>%
        group_by(sraIds, mut) %>%
        summarise(gtexIds_samples = gtexIds_samples[1], tissue = tissue[1], vaf_median = median(!!as.name(metric))) %>%
        ungroup() %>%
        spread(mut,vaf_median)
    
    x <- full_join(x, x_all, by='sraIds')
    
    return(x)
}

# Get median(VAF syn) / median(VAF non) for samples using eithermissense or nonsense as non
get_vaf_per_sample_impact <- function(x, metric) {
    
    vaf_per_sample_all <- x %>%
        group_by(sraIds,impact) %>%
        summarise(all=median(vaf)) %>%
        ungroup() %>%
        spread(impact, all) %>%
        mutate(all=Synonymous/!!as.name(metric)) %>%
        dplyr::select(sraIds,all)
    
    vaf_per_sample <- x %>%
        filter(!is.na(impact)) %>%
        group_by(sraIds, impact, mut) %>%
        summarise(gtexIds_samples = gtexIds_samples[1], tissue=tissue[1], vaf_median = median(vaf)) %>%
        ungroup() %>%
        group_by(mut) %>%
        spread(impact, vaf_median) %>%
        mutate(value = Synonymous/!!as.name(metric)) %>%
        filter(!is.na(value)) %>%
        ungroup() %>%
        dplyr::select(sraIds, gtexIds_samples, tissue, mut, value) %>%
        spread(mut, value)
    
    result <- full_join(vaf_per_sample, vaf_per_sample_all, by='sraIds')
    
    return(result)
}

# read covariates and do rankitNormalize
read_covariates <- function(metadata_file, sraIds, binary_covs='GENDER') {
    
    covariates <- read.table(metadata_file, header=T, sep="\t", stringsAsFactors=F) 
    covariates <- covariates[covariates$RACE == 3,]
    rownames(covariates)<-covariates$sraId
    covariates <- dplyr::select(covariates, -sraIds) 
    
    # Getting only for indicated sraIds
    covariates <- t(covariates[rownames(covariates) %in% sraIds,])
    covariates <- covariates[apply(covariates, 1, function(x) length(unique(x)) > 1),]
    
    covariates_a <- covariates[!rownames(covariates) %in% binary_covs,,drop=F]
    covariates_a <- rankitNormalize(covariates_a, IND=1)
    covariates_b <- covariates[rownames(covariates) %in% binary_covs,,drop=F]
    
    covariates <- rbind(covariates_a, covariates_b)
}


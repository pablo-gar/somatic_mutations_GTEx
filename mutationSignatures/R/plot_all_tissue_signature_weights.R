library("purrr")
library("tidyr")
library("dplyr")
library("ggplot2")
source("../../R/misc.R", chdir = T)

main <- function(cmdArgs = commandArgs(T)) {
    
    outfile <- cmdArgs[1]
    signature_orders <- as.numeric(parseArg(cmdArgs[2], sep = ","))
    H_consensus_files <- cmdArgs[-1]
    
    signature_orders <- as.numeric(parseArg("5,2,4,1,3", sep = ","))
    H_consensus_files <- list.files("/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/", full = T, pattern = "*H_cons*")
    
    H_all <- read_H_tables(H_consensus_files)
    
    if(any(!signature_orders %in% H_all$signature) | any(!H_all$signature %in% signature_orders))
        stop("Signture order labels do not match signatures")
    
    H_all$signature <- factor(H_all$signature, levels = signature_orders, ordered = T)
    ggplot(H_all, aes(x = signature, y = Weight)) +
        geom_jitter(aes(colour = signature), alpha = 0.3) +
        geom_boxplot(aes(fill = signature), width = 0.3) +
        facet_wrap(.~tissue, nrow = 5)+
        theme_bw()
    
}

read_H_tables <- function(x) {
    
    results <- map_dfr(x, function(y) {
                           tissue <- gsub("(.+?)-.*", "\\1", basename(y))
                           y <- read.table(y, sep = "\t", stringsAsFactors = F, header = T)
                           y <- y %>%
                               gather("Individual", "Weight", -signature) 
                           y$tissue <- tissue
                           return(y)
                          })
    
    return(results)
                               
    
}

main()

library("ggplot2")
source("../../R/ggthemes.R")

main <- function(cmdArgs = commandArgs(T)) {

    input_file <- cmdArgs[1]
    out_prefix <- cmdArgs[2]
    
    #input_file <- "/scratch/users/paedugar/somaticMutationsProject/supp_info/tables/Table_mutations_all.tsv.gzip"
    #out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/MAFhistograms_all/"

    mutations <- read.table(gzfile(input_file), sep = "\t", stringsAsFactors=F, header=T)
    
    mutations$VAF <- mutations$alt_count / mutations$coverage
    stats <- as.data.frame(quantile(mutations$VAF, probs = seq(0,1, length.out=5)))
    
    p <- ggplot(mutations, aes(x=VAF)) +
        geom_density() +
        scale_x_continuous(breaks=seq(0,0.7, by=0.1)) + 
        theme_noGrid()
    
    p2 <- ggplot(mutations, aes(x=VAF)) +
        geom_histogram(colour="grey20", fill="grey80", bins=80) +
        scale_x_continuous(breaks=seq(0,0.7, by=0.1)) + 
        theme_noGrid()
    
    ggsave(paste0(out_prefix, "MAF_density.pdf"), p, width=5, height=5)
    ggsave(paste0(out_prefix, "MAF_hist.pdf"), p2, width=5, height=5)
    
    write.table(stats, paste0(out_prefix, "stats.txt"), sep="\t", quote=F, col.names=F)
}

main()


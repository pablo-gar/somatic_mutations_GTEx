# Usage
# Rscript mutationVsMetadata_BiologicalEffects.R outBar.pdf readCountTable.txt n6_0.0_0.5 tissueCountFolder1 [tissueCountFolder1] [...] 
   

if(!dir.exists("../../R"))
   stop("Can't find the share R path, make sure to run script from directory where it lives")

library("ggplot2")
library("reshape")
library("purrr")
library("dplyr")
library("tidyr")
source("../../R/plots.R")
source("../../R/ggthemes.R")
source("../../R/misc.R")
source("../../R/gtex.R", chdir = T)

args <- commandArgs(T)
outPlotPrefix <- args[1]
outPlot_countsPrefix <- args[2]
outPrefix_corPlots <- args[3]
readCountFile <- args[4]
transDiversityFile <- args[5]
maf <- args[6]
folders <- args[7:length(args)]

#readCountFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/readCounts.txt_OLD"
#maf <- "n6_0.0_0.7"
#workingDir <- "/scratch/users/paedugar/somaticMutationsProject/clusterMutations/nonCluster/mutationCount/count/"
#workingDir <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/count/"
#transDiversityFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/transcriptome_shannon_diversity.txt"
#folders <- list.dirs(workingDir, recursive = F)
#outPrefix_corPlots <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutationsVsMetadata/2_spearman_correlations-n6_0.0_0.7"
#outPlotPrefix <- "~/out"
#outPlot_countsPrefix <- "~/out_counts"
folders <- file.path(folders, maf)
folders <- folders[!grepl("EXO", folders)]


# For bar plot
BAR_WIDTH <- 0.8
BAR_INCHES <- 0.7
PLOT_WIDTH <- 4

# For count plot
PLOT_HEIGHT_COUNT <- 12
PLOT_WIDTH_COUNT <- 4

# number of mutation types
N_MUT <- 7

metadataCols <- c("SMRIN", #RIN
                  "SMTSISCH", #Ischemic time
                  "SMTSPAX" #Time in fixative
                  )

indInfo <- c("AGE", "GENDER", "BMI")

factorVariables <- c("GENDER")
featureNames <- c("RIN", "Ischemic time", "Time in fixative", "Age", "Self-defined race", "Gender", "Uniquely mapped reads", "Transcriptome diversity", "Genotype PC1", "Genotype PC2", "Genotype PC3", "BMI")
names(featureNames) <- c("SMRIN", "SMTSISCH", "SMTSPAX", "AGE", "ETHNCTY", "GENDER", "n_uniqueMapped", "transcriptome_diversity", "C1", "C2", "C3", "BMI")

artifacts <- c("SMRIN", "SMTSISCH", "SMTSPAX",  "n_uniqueMapped", "transcriptome_diversity")
noArtifacts <- names(featureNames)[!names(featureNames) %in% artifacts]

#------------------------------------#
# METHODS

# Read count table should be produced with ../bin/getReadCounts
readGeneCounts <- function(x) {
    readCountTable <- read.table(x, sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
    return(readCountTable)
}

readTransDiversity <- function(x) { 
    transDiversity <- read.table(x, sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
    return(transDiversity)
}

# Read mutation count tables
readMutations <- function(files) {
    countTable <- data.frame(sample = "", mutations_all = rep(0, length(args)), T.A = 0, T.G = 0, T.C = 0, C.A = 0, C.T = 0, C.G = 0, stringsAsFactors = F)

    for (i in 1:length(files)) {
        mutation <- read.table(files[i], sep = "\t", stringsAsFactors = F, header = F)
        countTable[i, "sample"] <- gsub(".txt", "", basename(files[i]))
        countTable[i, "T.A"] <- mutation[1,1]
        countTable[i, "T.G"] <- mutation[1,3]
        countTable[i, "T.C"] <- mutation[1,4]
        countTable[i, "C.A"] <- mutation[2,1]
        countTable[i, "C.T"] <- mutation[2,2]
        countTable[i, "C.G"] <- mutation[2,3]
        countTable[i, "mutations_all"] <- sum(mutation)
    }
    
    return(countTable)
    
}

get_permutated_mean <- function(resids, noArtifactsTable, realCors, realPvals, oneSided = T, nPerm = 10000) {
    permResults <- list()
    permPvals <- list()
    realPvals[is.na(realCors)] <- 1
    realCors[is.na(realCors)] <- 0
    
    for(i in 1:nPerm) {
        randomCors <- cor(resids[sample(length(resids))], noArtifactsTable, method = "sp")
        permResults[[i]] <- realCors - randomCors
        #permPvals[[i]] <- abs(randomCors) >= abs(realCors) 
        if(oneSided) {
            permPvals[[i]] <- rep(F, length(realCors))
            if(any(realCors > 0))
                permPvals[[i]][realCors > 0] <- randomCors[realCors > 0] >= realCors[realCors > 0]
            if(any(realCors <= 0))
                permPvals[[i]][realCors <= 0] <- randomCors[realCors <= 0] <= realCors[realCors <= 0]
        } else {
            permPvals[[i]] <- abs(randomCors) >= abs(realCors) 
        }
    }
    
    permResults <- do.call(rbind, permResults)
    permPvals <- do.call(rbind, permPvals)
    
    
    means <- colMeans(permResults)
    sds <- apply(permResults, 2, sd)
    pvals <- colSums(permPvals) / nrow(permPvals)
    
    return(data.frame(measure = names(means), realCor = as.vector(realCors), realPval = as.vector(realPvals), mean = means, sd = sds, pval = pvals, stringsAsFactors = F))
    
    
        
}

plotCorrelations <- function(correlationTablePerm, outPrefix) {
    
    tissueOrders <- 
        correlationTablePerm %>%
        group_by(measure) %>%
        mutate(tissueOrd = tissue[order(realCor)], rank = 1:length(tissueOrd)) %>%
        select(tissueOrd, rank) %>%
        ungroup() %>%
        spread(measure, tissueOrd)
                                        
    for(m in unique(correlationTablePerm$measure)) {

        toPlot <- 
            correlationTablePerm %>%
            filter(!is.na(mean), measure == m) %>%
            mutate(tissue = factor(tissue, levels = tissueOrders[,m, drop = T], ordered = T)) %>%
            mutate(pval_fdr = p.adjust(pval, method = "BH")) %>%
            mutate(signif_label = case_when(pval_fdr < 0.2 & pval_fdr >= 0.1 ~ "*",
                                            pval_fdr < 0.1 & pval_fdr >= 0.05  ~ "**", 
                                            pval_fdr < 0.05 & pval_fdr >= 0.01  ~ "***", 
                                            pval_fdr < 0.01 ~ "****",
                                            TRUE ~ ""))

        #p <- ggplot(toPlot, aes(x = tissue, y = mean)) + 
                #geom_pointrange(aes(y = mean, ymin = mean - sd, ymax = mean + sd)) +
        p <- ggplot(toPlot, aes(x = tissue, y = realCor)) + 
                geom_point() +
                geom_hline(yintercept = 0, colour = "grey80",  linetype = "dashed") + 
                annotate("rect", ymax = 0.05, ymin = -0.05, xmax = Inf, xmin = -Inf, alpha = 0.1, fill = "blue" ) +
                geom_text(aes(x = tissue, y = realCor + 0.02, label = signif_label), hjust = 0) +
                #geom_text(aes(label = paste0("rho = ", signif(realCor,2))), vjust = 0, size = 3) +
                coord_flip(ylim = c(-(max(abs(toPlot$mean))) - 0.15, max(abs(toPlot$mean)) + 0.15)) + 
                theme_bw()
            
        outname <- paste0(outPrefix, m, ".pdf")
        ggsave(outname, p, width = 4.8, height = 5.5)

        outname <- paste0(outPrefix, m, "_correlation_table.txt")
        write.table(toPlot, outname, sep = "\t", col.names = T, row.names = F, quote = F)
    }
    
}
    


#------------------------------------#
# MAIN

correlationTables <- list()
correlationTablesPerm <- list()
countTableAll <- list()
allArtResid <- list()
pvalues <- list()
all_artResiduals <- list()
counter <- 1
for (tissueFolder in folders) {
    
    tissue <- basename(dirname(tissueFolder))
    cat("Working with tissue", tissue, "\n")
    
    args <- file.path(tissueFolder, list.files(tissueFolder))

    # Read mutation count tables
    countTable <- readMutations(args)
    sraIds <- countTable$sample
    gtexIds <- sraToGtex(countTable$sample, formatOut = "short")
    gtexIdsLong <- sraToGtex(countTable$sample, formatOut = "long")

    rownames(countTable) <-  gtexIdsLong

    # Read gene counts file
    readCountTable <- readGeneCounts(readCountFile)
    
    # Read transcriptome diversity
    transDiversity <- readTransDiversity(transDiversityFile)

    # Read metadata
    sampleMeta <- readSampleAnnotationGTEX(metadataCols)
    indMeta <- readMetadataGTEX(indInfo)
    genotypePCA <- t(readGenotypePCA())


    # Getting ids for which we have info for all
    allInfo <- sraIds %in% rownames(readCountTable) & 
        gtexIdsLong %in% rownames(sampleMeta) & 
        gtexIds %in% rownames(indMeta) & 
        gtexIds %in% rownames(genotypePCA) &
        gtexIdsLong %in% rownames(transDiversity)

    sraIds <- sraIds[allInfo]
    gtexIds<- gtexIds[allInfo]
    gtexIdsLong <- gtexIdsLong[allInfo]

    # Merging tables
    countTable <- countTable[gtexIdsLong,]

    countTable <- cbind(countTable, readCountTable[sraIds,"n_uniqueMapped", drop = F])
    countTable <- cbind(countTable, transDiversity[gtexIdsLong, , drop = F])
    countTable <- cbind(countTable, sampleMeta[gtexIdsLong,])
    countTable <- cbind(countTable, indMeta[gtexIds,])
    countTable <- cbind(countTable, genotypePCA[gtexIds,1:3])
    countTable <- countTable[,-1]

    # Eliminate columns with all na
    countTable <- countTable[,colSums(is.na(countTable)) != nrow(countTable)]

    # Doing lm and getting pvalues for each coefficient

    if(nrow(countTable) > 1) {
        
        for(currentMut in colnames(countTable)[1:N_MUT]) {
            
            countTable$mutations <- rankitNormalize_vector(countTable[,currentMut])
        
            # Getting only biological data
            # Regressing out artifacts
            cModel <- lm(mutations ~ ., data = countTable[,colnames(countTable) %in% c("mutations", artifacts)])
            artResiduals <- resid(cModel)
            countTable <- countTable[names(artResiduals),]
            countTable$artResiduals <- artResiduals
            
            # Getting prediction of mutation values
            m <- fitted(cModel)
            s <- summary(cModel)$sigma
            average_simulations <- do.call(cbind, lapply(1:1000, function(x) rnorm(length(m), m, s)))
            allArtResid [[counter]] <- data.frame(mutation_residual = rowMeans(average_simulations), tissue = tissue, mut_type = currentMut, stringsAsFactors = F)
            
            countTableAll[[counter]] <- countTable
            countTableAll[[counter]]$tissue <- tissue
            countTableAll[[counter]]$mut_type <- currentMut
            
            # Getting correlations
            resids <- countTable$artResiduals
            noArtifactsTable <- countTable[, colnames(countTable) %in% noArtifacts]
            originalCor <- cor(resids, noArtifactsTable, method = "sp")
            originalCor_pvals <- sapply(noArtifactsTable, function(x) cor.test(resids, x, method = "sp")$p.value)
            correlationTables[[counter]] <- as.data.frame(originalCor)
            correlationTables[[counter]]$tissue <- tissue
            correlationTables[[counter]]$mut_type <- currentMut
            
            # Permutation
            correlationTablesPerm[[counter]] <- get_permutated_mean(resids, noArtifactsTable, originalCor, originalCor_pvals, oneSided = T)
            correlationTablesPerm[[counter]]$tissue <- tissue
            correlationTablesPerm[[counter]]$mut_type <- currentMut
            
            # Getting pvalues
            pvaluesCurrent <- summary(lm(artResiduals~., data = countTable[, colnames(countTable) %in% c("artResiduals", noArtifacts)]))$coefficients[,"Pr(>|t|)"]
            pvaluesCurrent[noArtifacts[!noArtifacts %in% names(pvaluesCurrent)]] <- 1
            pvaluesCurrent <- pvaluesCurrent[noArtifacts]
            pvalues[[counter]] <- as.data.frame(t(pvaluesCurrent))
            pvalues[[counter]]$tissue <- tissue
            pvalues[[counter]]$mut_type <- currentMut
            
            counter <- counter + 1
        }
    }
    
}


# Gathering all art residuals
allArtResid <- do.call(rbind, allArtResid)

toNumeric <- function(x) as.numeric(unlist(strsplit(x, ",")))

medianArtResid <- 
    allArtResid %>%
    group_by(tissue, mut_type) %>%
    summarise(med = median(mutation_residual), confidence = bootstrap_confidence_interval(mutation_residual, median, 10000, interval = 0.9, out_string = T), low_ci = toNumeric(confidence)[1], high_ci = toNumeric(confidence)[2]) %>%
    #mutate(low_ci = unlist(strsplit(confidence, ","))[1], high_ci = unlist(strsplit(confidence, ","))[2])
    ungroup()
tissueOrds <- medianArtResid[medianArtResid$mut_type == "mutations_all",]
medianArtResid$tissue <- factor(medianArtResid$tissue, levels = tissueOrds[order(tissueOrds$med), "tissue", drop = T], ordered = T)

# Gather all correlations
correlationTable <- as.data.frame(do.call(rbind, correlationTables))
correlationTable <- melt.data.frame(correlationTable, id.vars = c("tissue", "mut_type"))
rownames(correlationTable) <- paste0(correlationTable$tissue, ".", correlationTable$variable, ".",  correlationTable$mut_type)
correlationTable$variable <- featureNames[as.character(correlationTable$variable)]

# Gather permutated correlations
correlationTablePerm <- as.data.frame(do.call(rbind, correlationTablesPerm))

pvalue <- as.data.frame(do.call(rbind, pvalues))
pvalue <- melt.data.frame(pvalue, id.vars = c("tissue", "mut_type"))
pvalue <- pvalue[ !is.na(pvalue$value),]
rownames(pvalue) <- paste0(pvalue$tissue, ".", pvalue$variable, ".", pvalue$mut_type)



for (current_mut_type in unique(correlationTable$mut_type)) {
    
    pvalue_current <- pvalue[pvalue$mut_type == current_mut_type,]
    correlationTable_current <- correlationTable[ correlationTable$mut_type == current_mut_type,]
    correlationTablePerm_current <- correlationTablePerm[ correlationTablePerm$mut_type == current_mut_type,]
    medianArtResid_current <- medianArtResid[ medianArtResid$mut_type == current_mut_type,]
    
    
    toPlot <- correlationTable_current[rownames(pvalue_current[pvalue_current$value < 0.05,]), ]
    
    if(nrow(toPlot) > 0 ) {

        # Plot bars of significant metadata
        p <- ggplot(toPlot, aes(y = value, x = variable)) +
            geom_bar(stat = "identity", width = 0.5) +
            coord_cartesian(ylim = c(-0.35, 0.35)) +
            coord_flip() +
            facet_grid(tissue~., scales = "free_y") +
            ylab("Spearman correlation with somatic mutations") +
            xlab("")+
            theme_grid_x()
        
        ROWS <- nrow(toPlot)
        
        
    } else {
        
        p <- ggplot(data.frame(x = 1, y = 1), aes(x = x, y = y)) + geom_point()
        ROWS <- 1
        
    }
    
    # Plot and save all correlation for all metadata
    plotCorrelations(correlationTablePerm_current, paste0(outPrefix_corPlots, ".", current_mut_type, "."))
    
    # Plots residual across tissues
    p_resid <- ggplot(medianArtResid_current, aes(x = tissue, y = med)) +
        geom_pointrange(aes(ymin = low_ci, ymax = high_ci)) +
        coord_flip() +
        theme_bw()
    

    ggsave(paste0(outPlotPrefix, ".", current_mut_type, ".pdf"), p, width = PLOT_WIDTH, height = BAR_INCHES * ROWS)
    ggsave(paste0(outPlot_countsPrefix, ".", current_mut_type, ".pdf"), p_resid, width = PLOT_HEIGHT_COUNT, height = PLOT_WIDTH_COUNT)
}

# takes a table containing at least the columns, SRR sample id, tissue,
# It changes the SS sampled id to a new one and adds a subject id
source("../../R/gtex.R", chdir = T)
source("../../R/misc.R", chdir = T)

main <- function(cmdArgs = commandArgs(T)) {
    
    input <- cmdArgs[1]
    srr_index <- as.numeric(cmdArgs[2])
    colNames <- parseArg(cmdArgs[3], sep = ",")
    output <- cmdArgs[4]
    
    #input <- "/scratch/users/paedugar/somaticMutationsProject/supp_info/tables/Table_mutations_all.tsv.gzip.temp"
    #colNames <- c("chr","pos","ref","alt","context","coverage","alt_count","tissue","sample_id","subject_id")
    #srr_index <- 9
    
    if(colNames != "NULL") {
        tab <- read.table(input, sep = "\t", stringsAsFactors = F, header = F)
    } else {
        tab <- read.table(input, sep = "\t", stringsAsFactors = F, header = T)
    }
    
    # Create conversion vector fro sraIds
    sraIds_all <- read.table(CONFIG$auxiliaryFiles$sraTable, sep = "\t", stringsAsFactors = F, header =T)[,16]
    sraIds_uniq <- unique(sraIds_all)
    sraIds <- setNames(add_lead_zeroes(1:length(sraIds_uniq)), sraIds_uniq)
    
    # Create conversion data.frame for subject ids
    getxIds_long <- sraToGtex(sraIds_uniq, formatOut = "long")
    getxIds <- gsub("(GTEX-.+?)-.+", "\\1", getxIds_long)
    uniq_gtexIds <- unique(getxIds)
    uniq_gtexIds <- setNames(add_lead_zeroes(1:length(uniq_gtexIds)), uniq_gtexIds)
    gtexIds <- data.frame(sraIds = names(getxIds), gtexIds = getxIds, gtexIds_samples = getxIds_long, gtexIds_de = uniq_gtexIds[getxIds], stringsAsFactors = F)
    rownames(gtexIds) <- gtexIds$sraIds
    
    # Convert
    sraIds_original <- tab[,srr_index]
    tab[,srr_index] <- sraIds[sraIds_original]
    tab$subject <- gtexIds[sraIds_original,"gtexIds_de"]
    
    
    if(colNames != "NULL")
        colnames(tab) <- colNames
    
    tab_all <- tab
    tab_all <- cbind(tab_all, gtexIds[sraIds_original,])
    
    write.table(tab, output, sep = "\t", quote = F, col.names = T, row.names = F)
    write.table(tab_all, paste0(output, ".original"), sep = "\t", quote = F, col.names = T, row.names = F)
    
    
}

add_lead_zeroes <- function(x) {
    sprintf(paste0("%0", max(nchar(x)), "d"), x)
}

main()

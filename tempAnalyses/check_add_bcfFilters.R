newPileups = "/scratch/users/paedugar/somaticMutationsProject/pileups"
newPileup_base = "pileups"
oldPileup_base = "pileupsNoFilterStats"


allFiles <- file.path(newPileups, list.files(newPileups, recursive = T))
allFiles <- allFiles[!grepl("Muscle_Skeletal", allFiles)]
allFiles <- allFiles[!grepl("Spleen", allFiles)]

new_size <- rep(0, length(allFiles))
old_size <- new_size
for(i in seq_along(allFiles)) {
    
    this <- allFiles[i]
    old <- gsub(newPileup_base, oldPileup_base, this)
    
    flush.console()
    cat(this, "\n")
    
    if(new_size[i] != 0)
        next
    
    new_size[i] <- as.numeric(system(paste("wc -l <", this), intern = T))
    old_size[i] <- as.numeric(system(paste("wc -l <", old), intern = T))
}

boxFiles <- gsub(newPileups, "mappedReads", allFiles[new_size != old_size])
boxFiles <- gsub(".txt", "_RmdupSortedAligned.out.bam", boxFiles)
foldersToCreate <- unique(dirname(boxFiles))
boxFiles <- paste(boxFiles, "-o", boxFiles)
boxFiles <- paste(boxFiles, collapse = " ")

foldersToCreate <- paste("mkdir -p", paste(foldersToCreate, collapse = " "))
writeLines(boxFiles, "boxFiles")
writeLines(foldersToCreate, "foldersToCreate")

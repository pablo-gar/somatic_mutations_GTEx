library("purrr")

map_files <- list.files("/scratch/users/paedugar/somaticMutationsProject/mutationCount/filtered/Nerve_Tibial/n6_0.0_0.7", full = T)[1:2]

map_files <- "mutations_wihtGermLine_Nerve.txt"
map_files <- "mutations_wihtGermLine_Nerve_n1_c10_0.0_1.0.txt"
map_files <- "mutations_wihtGermLine_Nerve_n6_c40_0.0_1.0.txt"
map_files <- "mutations_wihtGermLine_Nerve_n4_c20_0.0_1.0.txt"

map_mut <- map_dfr(map_files, function(x) read.table(x, sep = "\t", stringsAsFactors = F, header = F))
map_mut$VAF <- map_mut$V7 / map_mut$V6

hist(map_mut$VAF, breaks = 200)
# The dip is already there

# Look to  see if there is anything special for mutations close to VAF 0.5
map_mut$vaf_group <- "rest"

# Great groups of vars
VAFS <-  seq(0.2, 0.7, by = 0.05)
for(i in VAFS) {
    map_mut$vaf_group[map_mut$VAF > (i - 0.02) & map_mut$VAF < (i + 0.02)] <- as.character(i)
}

# Anything special about n support?
toPlot <- map_mut[map_mut$vaf_group != "rest",]
plot(x = VAFS, y = tapply(toPlot$V7, toPlot$vaf_group, mean), ylab = "average reads variant", pch = 19)
plot(x = VAFS, y = tapply(toPlot$V6, toPlot$vaf_group, mean), ylab = "average coverage", pch = 19)




#######

# Plot density after eliminating each filter
filters <- unique(unlist(strsplit(map_mut$V8, ";")))
filters <- filters[!filters %in% c("PASS", "clustered_mutation")]
for(i in filters) {
    
    dev.new()
    plot(density(map_mut[grep(i, map_mut$V8), "VAF"]), main = i)
}

plot(density(map_mut[ map_mut$V8 %in% c("PASS", "clustered_mutation"), "VAF"]))
    

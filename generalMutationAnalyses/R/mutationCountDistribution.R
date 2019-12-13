library(ggplot2)

args <- commandArgs(T)

outplot <- args[1]
args <- args[-1]

countTable <- data.frame(sample = "", mutations = rep(0, length(args)), stringsAsFactors = F)
for (i in 1:length(args)) {
    mutation <- read.table(args[i], sep = "\t", stringsAsFactors = F, header = F)
    countTable[i, "sample"] <- basename(args[i])
    countTable[i, "mutations"] <- sum(mutation)
}


stats <- data.frame(x = Inf, y = Inf,
                    label = paste0(
                        "nInds = ", nrow(countTable), "\n",
                        "mean = ", mean(countTable$mutations), "\n",
                        "sd = ", sd(countTable$mutations), "\n",
                        "min = ", min(countTable$mutations), "\n",
                        "max = ", max(countTable$mutations), "\n"
                        )
                        
                    )
                    

p <- ggplot(countTable, aes( x = mutations)) +
geom_histogram(bins = floor(nrow(countTable) / (nrow(countTable) * 0.05))) +
geom_text(aes(x = x, y = y, label = label), data = stats, vjust = 1, hjust = 1) +
theme_bw()

ggsave(outplot, p)

# Plots the theoritcal probabilities of observing a mutation
# as an artifact of sequencing errors giving the parameters of
# sequencing error, coverage and n_suport

library("ggplot2")
library("dplyr")


cmdArgs <- commandArgs(T)
outFile <- cmdArgs[1]

r_cutoff <- 8
exome_coverage <- c(20, 40, 80, 160)
coverage <- c(40, 100, 200, 400, 800, 1600, 2400, 4800)
n <- 6:39
#qs <- 30:40
#base_multiplier <- 1/3 # this is a multiplier to the probability of qs that accounts for having three possible bases to mutate to
#qs_convert <- function(x) 10^(-x/10)

probs <- expand.grid(coverage, n, exome_coverage)
colnames(probs) <- c("coverage", "n", "exome_coverage")
probs$r <- probs$n / probs$coverage * probs$exome_coverage

#minProb <- 
#    probs %>%
#    group_by(coverage) %>%
#    summarise(label = paste0("max(p) = ", signif(max(prob), 2)), n = -Inf, prob = Inf, qs = min(qs)) %>%
#    ungroup()

p <- ggplot(probs, aes(y=r, x = n)) +
    geom_line(aes(colour = as.factor(coverage), group = coverage)) +
    scale_colour_manual(values = rainbow(length(coverage), s = 0.5)) +
    geom_hline(yintercept = r_cutoff, linetype = "dashed", colour = "grey40") +
    labs(colour = "Coverage RNA-seq") + 
    #geom_text(aes(label = label), data = minProb, vjust = 1, hjust = 0) + 
    scale_y_continuous(trans = "log10") +
    xlab("Reads supporting ALT allele") + 
    ylab("r") + 
    facet_wrap(exome_coverage~., scales = "free", ncol = 2)+
    theme_bw()

ggsave(outFile, p, width = 7, height = 7)



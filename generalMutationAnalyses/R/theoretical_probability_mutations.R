# Plots the theoritcal probabilities of observing a mutation
# as an artifact of sequencing errors giving the parameters of
# sequencing error, coverage and n_suport

library("ggplot2")
library("dplyr")


cmdArgs <- commandArgs(T)
outFile <- cmdArgs[1]

#coverage <- round(seq(40, 8000, length.out = 100))
coverage <- c(40, 10e3)
n <- 6:39
qs <- 30:40
base_multiplier <- 1/3 # this is a multiplier to the probability of qs that accounts for having three possible bases to mutate to
qs_convert <- function(x) 10^(-x/10)

probs <- expand.grid(coverage, n, qs)
colnames(probs) <- c("coverage", "n", "qs")
probs$prob <- pbinom(probs$n, probs$coverage, qs_convert(probs$qs) * base_multiplier, lower.tail = F)

minProb <- 
    probs %>%
    group_by(coverage) %>%
    summarise(label = paste0("max(p) = ", signif(max(prob), 2)), n = -Inf, prob = Inf, qs = min(qs)) %>%
    ungroup()

p <- ggplot(probs, aes(y=prob, x = n)) +
    geom_line(aes(colour = as.factor(qs), group = qs)) +
    scale_colour_manual(values = heat.colors(length(qs))) +
    scale_y_continuous(trans = "log10") +
    labs(colour = "Seq qual score") + 
    xlab("Reads supporting ALT allele") + 
    geom_text(aes(label = label), data = minProb, vjust = 1, hjust = 0) + 
    ylab("p(sequencing error)") + 
    facet_wrap(coverage~., scales = "free", ncol = 1)+
    theme_bw()

ggsave(outFile, p, width = 7, height = 7)



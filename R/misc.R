library(dplyr)

printTime <- function(x = "", carriageReturn = F) {
	# Prints to console the current date and time followed by
	# the message(s) in x
	
	first = ""
	if(carriageReturn) first = "\r"
	
	flush.console()
	cat(first, as.character(Sys.time()), x)
}

#' Calcualtes bootstrap confidence interval after applying a function to a vector
#' @param x vector
#' @param FUN function to be applied to vector, has to return a numeric vector of size 1
#' @param bootstrap_counts number of bootrsap rounds
#' @param interval value between 0 and 1 depecting the confidence interval desired
#' @param bootstrap_size value between 0 and 1 indicanting the percentage of data taken at each bootstrap repetition
#' @param out_string logical, if true returns a string representation of the interval, otherwise return a vector of size 2
bootstrap_confidence_interval <- function(x, FUN, bootstrap_counts = 1000, interval = 0.95, bootstrap_size = 1, out_string = F) {
    
    bootstrap_counts <- as.integer(bootstrap_counts[1])
    interval <- interval[1]
    bootstrap_size <- bootstrap_size[1]
    out_string <- out_string[1]
    
    if(!is.numeric(x) | !is.vector(x))
        stop("x has to be a numeric vector")
    
    if(bootstrap_counts < 1)
        stop("bootstrap_counts has to be a positive integer")
    
    if(interval < 0 | interval > 1)
        stop("interval has to be between 0 and 1")
    
    if(bootstrap_size < 0 | bootstrap_size > 1)
        stop("bootstrap_size has to be between 0 and 1")
    
    if(!is.logical(out_string))
        stop("out_string has to be logical")
    
    if(!is.function(FUN))
        stop("FUN has to be a function")
    
    original <- FUN(x)
    
    results <- rep(NA, bootstrap_counts)
    bSize <- round(length(x) * bootstrap_size)
    for(i in 1:bootstrap_counts)
        results[i] <- original - FUN(x[sample(bSize, replace = T)])
    
    qSize <- (1 - interval) / 2 
    confidence_interval <- original + quantile(results, c(0 + qSize, 1 - qSize))
    
    if (out_string)
        confidence_interval <- paste(confidence_interval, collapse = ",")
    
    return(confidence_interval)
    
}

#' Gets labels of signifance for a vector of pvalues
#' @param pvalues vector of pvalues
#' @param label vector of the labels corresponding to each interval in `cutoff`
#' @param cutoff vector of increasing cutoffs
#' @param na.rm logical, if TRUE makes all NAs in the `pvalues` vectors into 1s
labelPvalues <- function(pvalues, label = c("***", "**", "*", "n.s."), cutoff = c(0, 0.001, 0.01, 0.05, 1), na.rm = T) {
    
    
    if(length(cutoff) - 1 != length(label))
        stop("the label vector has to be 1 - length of cutoff vector")
    
    if(na.rm)
        pvalues[is.na(pvalues)] <- 1
    
    if(any(pvalues < 0 | pvalues > 1))
        stop("Pvalues have to be between 0 and 1")
    
    pvalLabels <- label[findInterval(pvalues, cutoff)]
    pvalLabels[is.na(pvalLabels)] <- label[length(label)]
    
    return(pvalLabels)
}

parseArg <- function(x, sep, trim = T) {
		
	#Takes a string an creates a vector using sep as the separator
	# x - a string, if a vector only first element will be considred
	# sep - a string, if a vector only first element will be considered
	# trim - bool, if true white spaces will be eliminated
		
	if(!is.character(x))
		stop ("x has to be a string")
	if(!is.character(sep))
		stop ("x has to be a string")
	
	x <- x[1]
	sep <- sep[1]
	
	if(trim)
		x <- trim(x)
	
	return (unlist(strsplit(x,sep)))
	
}

trim <- function(x) { 
		
	# Trims white space from a character vector
	
	if(!is.character(x))
		stop("x has to be string or a character vector")

	return(gsub("\\s+", "", x))
	
}

grepTempFile <- function(x, pattern, tempLocation = "."){
	        
	# Creates a new file based of x only with lines containing
	# the specified pattern(s)
	#
	# WARNING: uses unix grep
	#   
	# x - string - path to file
	# pattern - vector - patterns to select in file

	patternArg <- paste0("-e ", paste0(pattern, collapse = " -e "))
	outFile <- file.path(tempLocation, paste0(basename(x), ".temp", sample(1:1000, 1)))

	system(paste("grep", patternArg, x, ">", outFile))

	return(outFile)

}

rankitNormalize <- function(x, IND = 1) {

    # Normalizes rows (IND = 1) or columns (IND = 2)
    # to a quantile standard normalization (i.e. rankit)

    stopifnot(is.matrix(x))
    stopifnot(is.numeric(x))
    
    rowNames <- rownames(x)
    colNames <- colnames(x)

    x <- apply(x, IND, rankitNormalize_vector)
    if(IND == 1)
        x <- t(x)
    
    rownames(x) <- rowNames
    colnames(x) <- colNames

    return(x)

}

rankitNormalize_vector <- function(x) {

    stopifnot(is.numeric(x))
    noNa <- !is.na(x)
    x[noNa] <- qnorm((rank(x[noNa]) - 0.5) / sum(noNa))
    return(x)

}

sourceDir <- function (path, pattern = "\\.[rR]$", env = NULL, chdir = TRUE) {
    
    # Sources all files in a folder using relative paths 
    
    files <- sort(dir(path, pattern, full.names = TRUE))
    lapply(files, source, chdir = chdir)
}

#' Gets the slope of two points
#' @param x1 x coordinate of first point
#' @param x2 x coordinate of second point
#' @param y1 y coordinate of first point
#' @param y2 y coordinate of second point
slope <- function(x1, x2, y1, y2) {
    return( (y2-y1) / (x2-x1) )
}

#' Gets the y intercept of a point and slople
#' @param x x coordinate point
#' @param y y coordinate point
#' @param m slope 
y_intercept <- function(x, y, m) {
    return( y - m*x)
}

#' Creates a blue-black-yellow color vector intended for signed pvalues in log scale.
#' The size of the resulting vector is 1 + 2(_maxPval_ - significant + 1)
#' @param maxPval the max pvalue available
#' @param significant the lowest significant pvalue
get_colors_significant <- function(maxPval = 10, significant = 2) {

    blues <- colorRampPalette(c("#4949FF", "#2F2F50"))(maxPval - significant + 1)
    yellows <- colorRampPalette(c("#50502F", "#FFFF4A"))(maxPval - significant + 1)
    blacks <- rep("grey20", significant * 2 -1)

    return(c(blues,blacks,yellows))
}

#' Fisher's method to combine pvalues
#' @param x a vector of pvalues
#' @return a numeric pvalue

fisher.pvalues  <- function (x){
    
    stopifnot(is.numeric(x))
    stopifnot(all(x >= 0 & x <= 1))
    
    # To account for zeroes
    x <- x + 1e-15
    
    deg_free <- length(x) * 2
    y <- -2 * sum(log(x))
    results <- 1 - pchisq(y, df = deg_free);
    results <- as.numeric(results);
    return(results) 
}

rolling_median <- function(x, y = NULL, windows = 100, by_quantile = T) {
    
    if(is.null(y))
        y <- x
    
    if(by_quantile) {
        groups <- quantile(y, seq(0,1, length.out=windows))
    } else {
        groups <- seq(min(y), max(y), length.out = windows)
    }
    
    dat <- data.frame(x = x, y = y)
    
    dat$group <- findInterval(dat$y, groups, all.inside = T)
    
    medians <-
        dat %>%
        group_by(group) %>%
        summarise(median = median(x)) %>%
        ungroup() %>%
        as.data.frame()
    rownames(medians) <- medians$group
    
    dat$median_mut_group <- medians[as.character(dat$group), "median"]
    
    result <- dat$x - dat$median_mut_group
    
    return(result)
    
}

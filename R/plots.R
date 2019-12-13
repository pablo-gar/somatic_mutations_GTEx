library(ggplot2)
library(reshape)
library(gplots)

Heatmap <- function(x, heatmap_gradient = NULL, dendrogram = "none", Rowv = FALSE, Colv = FALSE, symbreaks = F, trace = "none", density = "none", ...) {
    
    library(gplots)
    
    if(is.null(heatmap_gradient)) 
        heatmap_gradient <- colorRampPalette(c("#FFF5F0", "#3399ff"))(20)
    
    x <- as.matrix(x)
    
    #----
    # Gathering pars for heatmap2
    heatPars <- list(x = x,
                     col = heatmap_gradient, 
                     dendrogram = dendrogram, Rowv = Rowv, Colv = Colv,
                     symbreaks = symbreaks,
                     #scale = "column",
                     trace = trace, density = density,
                     ...
                     )
    
    #----
    
    return(do.call(heatmap.2, heatPars))
    
    
}

pointRange <- function(dataframe, x, y, errorBarYmax = NULL, errorBarYmin = NULL,scales = "free",
                       
                       # Category by color
                       colour = NULL,
                       
                       # Label settings
                       labelSize = 4, labelRound = 2,
                       
                       # Facet settings
                       facet_x = NULL, facet_y = NULL, nrowFactor = 1, ncolFactor = 1,
                        
                       # Point settings
                       alpha = 1, pColour = "black", pSize = 2,
                       
                       # Axis setting
                       ylab = NULL, xlab = NULL
                       
                       ) {
    #--------------------------
    # Checking parameters 
    #--------------------------
    
    if(!is.data.frame(dataframe)) stop ("dataframe has to be a data.frame")
    if(!is.character(x) | !is.character(y) | length(x) > 1 | length(y) > 1) stop ("x and y have to be character vectors of length one")
    
    if(!isCol(x, dataframe) | !isCol(y,dataframe)) stop ("x and y have to be columns of dataframe")
    
    if(!is.null(colour))
       if(!isCol(colour, dataframe)) 
           stop("colour has to be a column of dataframe or null")
    
    if(!is.null(facet_x))
       if(!isCol(facet_x, dataframe)) 
           stop("facet_x has to be a column of dataframe or null")
    if(!is.null(facet_y)) 
        if(!isCol(facet_y, dataframe)) 
            stop("facet_y has to be a column of dataframe or null")
    
    if(!is.null(errorBarYmax))
       if(!isCol(errorBarYmax, dataframe)) 
           stop("errorBarYmax has to be a column of dataframe or null")
    
    if(!is.null(errorBarYmin)) 
        if(!isCol(errorBarYmin, dataframe)) 
            stop("errorBarYmin has to be a column of dataframe or null")
    
    if((is.null(errorBarYmax) & !is.null(errorBarYmin) ) | (is.null(errorBarYmin) & !is.null(errorBarYmax) ) )
        stop("errorBarYmax and errorBarYmin have to be both specified or NULL")
    
    #-----------------------------
    # Assemblying Facets
    #-----------------------------
    
    if (is.null(facet_y) & is.null(facet_x)) {
        facetForm <- ""
    } else if(is.null(facet_y)) {
        facetForm <- paste0("~", facet_x)
        ncol <- length(unique(dataframe[,facet_x]))
        nrow <- 1 * nrowFactor
    } else if(is.null(facet_x)) {
        facetForm <- paste0(facet_y, "~")
        ncol <- 1 * ncolFactor
        nrow <- length(unique(dataframe[,facet_y]))
    } else {
        facetForm <- paste0(facet_y, "~", facet_x)
        ncol <- length(unique(dataframe[,facet_x]))
        nrow <- length(unique(dataframe[,facet_y]))
    }
    
    #----------------------------
    # Making plots
    #----------------------------
    
    # Plot
    pointRangeArgs <- list()
    aesPointRange <- list()
    
    # Base plot
    if(is.null(errorBarYmax)) {
        p <- ggplot(dataframe, aes_string( x = x, y = y))
        aesPointRange <- c(aesPointRange, list(ymin = y, ymax = y))
    } else {
        p <- ggplot(dataframe, aes_string( x = x, y = y))
        aesPointRange <- c(aesPointRange, list( ymin = errorBarYmin, ymax = errorBarYmax))
    }
       
    # Colour and size
    if(!is.null(colour)) {
        aesPointRange <- c(aesPointRange, list(colour = colour))
        pointRangeArgs <- c(pointRangeArgs, list(position = position_dodge(width = 0.5)))
    } else {
        pointRangeArgs <- c(pointRangeArgs, list(colour = pColour))
    }
    
    # Gathering arguments for point range
    aesPointRange <- do.call(aes_string, aesPointRange)
    pointRangeArgs <- c(pointRangeArgs, list(mapping = aesPointRange, alpha = alpha, size = pSize))
    p <- p + do.call(geom_pointrange, pointRangeArgs)
    
    # Axes
    if(!is.null(ylab))
        p <- p + ylab(ylab)
    if(!is.null(xlab))
        p <- p + xlab(xlab)
    
    # Facets 
    if(!is.null(facet_x) | !is.null(facet_y)){
        facetForm <- as.formula(facetForm)
        p <- p + facet_wrap(facetForm, scales = scales, nrow = nrow)
    }
    

    #DONE
    return(p) 
        
    
    
}

#' Makes a scatter plot from a data frame with the following columns
#' x = string; column with data for x axis
#' y = string, column with data for y axis
#' facet_x = factor, column with factor to subdivide the data and plot in different columns
#' facet_y = factor, column with factor to subdivide the data and plot in different rows
#'
#' dataframe = data.frame, with the columns indicated above
#'
#' ... further parameters for facet_wrap
#'
#' returns a plot object of ggplot2
#' @examples
#' np <- 200 
#' testData <- data.frame(x = rnorm(np), category_A = sample(c("A", "B"), np, replace = T), category_B = sample(c("X", "Y", "Z"), np, replace = T))
#' testData$y <- jitter(testData$x, amount = 1.5)
#' scatter(testData, x = "x", y = "y", facet_x = "category_A", facet_y = "category_B", regression = T)
scatter <- function(dataframe, x, y, scales = "free", 
                    
                    # Label settings
                    labelSize = 4, labelRound = 2, method_cor = "spearman",
                    
                    # Facet settings
                    facet_x = NULL, facet_y = NULL, nrowFactor = 1, ncolFactor = 1,
                     
                    # Point settings
                    alpha = 0.75, pColour = "black", pSize = 2,
                    
                    # Axis setting
                    ylab = NULL, xlab = NULL,
                    
                    # Regression line options
                    regression = F, regressionColour = "lightskyblue3", regressionSE = F
                    ) {
    
    
    #--------------------------
    # Checking parameters 
    #--------------------------
    
    if(!is.data.frame(dataframe)) stop ("dataframe has to be a data.frame")
    if(!is.character(x) | !is.character(y) | length(x) > 1 | length(y) > 1) stop ("x and y have to be character vectors of length one")
    if(!isCol(x, dataframe) | !isCol(y,dataframe)) stop ("x and y have to be columns of dataframe")
    
    if(!is.null(facet_x))
       if(!isCol(facet_x, dataframe)) 
           stop("facet_x has to be a column of dataframe or null")
    if(!is.null(facet_y)) 
        if(!isCol(facet_y, dataframe)) 
            stop("facet_y has to be a column of dataframe or null")
    
    
    #-----------------------------
    # Getting correlation strings
    #-----------------------------
    
    Params <- list(dataframe = dataframe, x = x, y = y, labelRound = labelRound, method = method_cor)
    
    if (is.null(facet_y) & is.null(facet_x)) {
        Params <- c(Params, list(category_col = NULL))
        facetForm <- ""
    } else if(is.null(facet_y)) {
        Params <- c(Params, list(category_col = facet_x))
        facetForm <- paste0("~", facet_x)
        nrow <- 1 * nrowFactor
        ncol <- ceiling(length(unique(dataframe[,facet_x])) / nrow)
    } else if(is.null(facet_x)) {
        Params <- c(Params, list(category_col = facet_y))
        facetForm <- paste0(facet_y, "~.")
        ncol <- 1 * ncolFactor
        nrow <- ceiling (length(unique(dataframe[,facet_y])) / ncol)
    } else {
        Params <- c(Params, list(category_col = c(facet_x, facet_y)))
        facetForm <- paste0(facet_y, "~", facet_x)
        ncol <- length(unique(dataframe[,facet_x]))
        nrow <- length(unique(dataframe[,facet_y]))
    }
    
    corString <- do.call(getCorString, Params)
    
    corString$x <- -Inf
    corString$y <- Inf
    colnames(corString)[ (ncol(corString) - 1) : ncol(corString)] <- c(x,y)
    
    #----------------------------
    # Making plots
    #----------------------------
    
    p <- ggplot(dataframe, aes_string( x = x, y = y )) 
    
    
    if(!is.null(ylab))
        p <- p + ylab(ylab)
    if(!is.null(xlab))
        p <- p + xlab(xlab)
    
    if(pColour %in%  colnames(dataframe)) {
        p <- p + geom_point(aes_string(colour = pColour), alpha = alpha, size = pSize)
    } else {
        p <- p + geom_point(alpha = alpha, size = pSize, colour = pColour)
    }
    
    p <- p + geom_text(aes(label = text), data = corString, hjust = 0, vjust = 1, size = labelSize, fontface = "italic") + 
    theme_bw()
    
    if(regression)
        p <- p + geom_smooth(method = "lm", colour = regressionColour, se = regressionSE)
    
    if(!is.null(facet_x) | !is.null(facet_y)){
        facetForm <- as.formula(facetForm)
        p <- p + facet_wrap(facetForm, scales = scales, nrow = nrow)
    }
    

    #DONE
    return(p) 
        
}

isCol <- function(x, mat) {
        return (x %in% colnames(mat))
}

getCorString <- function(dataframe, x, y, category_col = NULL, labelRound = 2, method = "pearson") {
        
        # Calculates pearson correlation coefficient and pvalue between two columns
        # in a dataframe
        
        #--------------------------
        # Checking parameters 
        #--------------------------
        if(!is.character(x) | !is.character(y) | length(x) > 1 | length(y) > 1) stop ("x and y have to be character vectors of length one")
        if(!isCol(x, dataframe) | !isCol(y,dataframe)) stop ("x and y have to be columns of dataframe")
        
        if(!is.null(category_col)) {
            for(i in category_col) { 
                if(!isCol(i, dataframe)){ stop ("catgories have to be columns of dataaframe")}
            }
        } 
        #-------------------------
        # Getting correlations
        #-------------------------
        
        # Construct category list of factors
        if(!is.null(category_col)) {
            catList <- list()
            for(i in category_col){
                catList[[i]] <- dataframe[,i]
            }
        } else {
            catList <- rep(T, nrow(dataframe))
        }
        
        # Loop by categorie 
        results <- by(dataframe, catList, function(x, xVal, yVal, category_col) { 
                      # Get correlations
                      corR <- cor.test(x[,xVal], x[,yVal], method = method)
                      corText <- paste0("r = ", signif(corR$estimate,labelRound) , "\np = ", signif(corR$p.value, labelRound))
                      result <- data.frame(r = corR$estimate, p =  corR$p.value, text = corText, stringsAsFactors = F)
                      # Append category names
                      if(!is.null(category_col)) {
                          for(i in category_col){
                              result <- cbind (result, x[1,i], stringsAsFactors = F)
                          }
                          colnames(result)[ (ncol(result) - length(category_col) + 1) : ncol(result) ] <- category_col
                      }
                      return(result)
                   }, xVal = x, yVal = y, category_col = category_col)
        
        results <- do.call(rbind, results)
        return(results)
        
}


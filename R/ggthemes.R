library("ggplot2")
library("proto")

plot_example <- function() {
    
    ggplot(diamonds, aes(x = carat, y = price)) +
        geom_point(aes(colour = color)) +
        labs(title = "Diamond price depends on its carat", 
             subtitle = "Divided by cut",
             caption = "Example plot") +
        facet_wrap(~cut)
    
}

plot_example_2 <- function() {
    
    ggplot(diamonds, aes(x = carat, y = price)) +
        geom_point(aes(colour = color)) +
        labs(title = "Diamond price depends on its carat", 
             subtitle = "Divided by cut",
             caption = "Example plot")
    
}


theme_sleek <- function() {
    theme_bw() +
        theme(
              text = element_text(family = "sans", colour = "grey35"),
              
              panel.border  = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              
              axis.title.x = element_text(size = rel(1.05), face = "bold"),
              axis.title.y = element_text(size = rel(1.05), face = "bold"),
              axis.line = element_line(colour = "grey30", lineend = "square", size = 0.4),
              axis.ticks = element_line(colour = "grey35", lineend = "square", size = 0.3),
              
              legend.box.background = element_rect(fill = "transparent", colour = "grey30"),
              legend.background = element_blank(),#element_rect(fill = "transparent", colour = "grey30"),
              legend.key = element_blank(),
              legend.title = element_text(face = "bold"),
              
              panel.background = element_blank(),
              
              plot.background = element_blank(),
              plot.title = element_text(size = rel(1.2), colour = "gray25"),
              plot.subtitle = element_text(size = rel(1), colour = "gray25", face = "italic"),
              plot.caption = element_text(size = rel(0.8), colour = "gray25", hjust = 1),
              
              #strip.background = element_rect(fill = "white", colour = "white", linetype = 1, size = 0.1),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold.italic", colour = "grey55", size = 10),
              
              )
}

theme_sleek_grey <- function() {
    theme_sleek() +
        theme(
              plot.background = element_rect(fill = "grey90", colour = "grey30", size = rel (1.2)),
              panel.background = element_rect(fill = "grey85", colour = "transparent"),
              legend.box.background = element_rect(fill = "grey85", colour = "transparent"),
              strip.background = element_rect(fill = "grey82", colour = "transparent"),
              strip.text = element_text(face = "bold.italic", colour = "grey45", size = 10),
              axis.line = element_blank()
              )
}

theme_grid <- function() {
    theme(panel.grid.major = element_line(size = rel(0.5), colour = "grey75"))
}

theme_grid_y <- function() {
    theme(panel.grid.major.y = element_line(size = rel(0.5), colour = "grey75"))
}

theme_grid_x <- function() {
    theme(panel.grid.major.x = element_line(size = rel(0.5), colour = "grey75"))
}

theme_base <- function(base_size = 12) {
  structure(list(
    axis.line =         element_blank(),
    axis.text.x =       element_text(colour = NA,size = base_size * 0.8 , lineheight = 0.9, vjust = 1),
    axis.text.y =       element_text(colour = NA,size = base_size * 0.8, lineheight = 0.9, hjust = 1),
    axis.ticks =        element_line(colour = NA, size = 0.2),
    axis.title.x =      element_text(colour = NA,size = base_size, vjust = 1),
    axis.title.y =      element_text(colour = NA,size = base_size, angle = 90, vjust = 0.5),
    axis.ticks.length = unit(0.3, "lines"),
    axis.ticks.margin = unit(0.5, "lines"),
 
    legend.background = element_rect(colour=NA), 
    legend.key =        element_rect(colour = NA, ),
    legend.key.size =   unit(1.2, "lines"),
    legend.text =       element_text(colour = NA,size = base_size * 0.8),
    legend.title =      element_text(colour = NA,size = base_size * 0.8, face = "bold", hjust = 0),
    legend.position =   "right",
 
    panel.background =  element_rect(fill = NA, colour = NA), 
    panel.border =      element_rect(fill = NA, colour=NA), 
    panel.grid.major =  element_line(colour = "grey90", size = 0.2),
    panel.grid.minor =  element_line(colour = NA, size = 0.5),
    panel.margin =      unit(0.25, "lines"),
 
    strip.background =  element_rect(fill = NA, colour = NA), 
    strip.text.x =      element_text(colour = NA,size = base_size * 0.8),
    strip.text.y =      element_text(colour = NA,size = base_size * 0.8, angle = -90),
 
    plot.background =   element_rect(colour = NA),
    plot.title =        element_text(colour = NA,size = base_size * 1.2),
    plot.margin =       unit(c(1, 1, 0.5, 0.5), "lines")
  ), class = "options")
}


theme_noGrid <- function (base_size = 12, base_family = "") 
{
    theme_grey(base_size = base_size, base_family = base_family) %+replace% 
        theme(
                axis.line = element_line(colour = "grey60"),     
                axis.text = element_text(size = rel(0.8)), 
                axis.ticks = element_line(colour = "black"), 
                legend.key = element_rect(colour = "grey80"), 
                panel.background = element_rect(fill = "white", colour = NA), 
                panel.border = element_rect(fill = NA,colour = NA),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.major.x = element_blank(),                
                panel.grid.minor = element_blank(), 
                strip.background = element_rect(fill = "white", colour = "grey50", size = 0.5)
            )
}


theme_grid_y <- function (base_size = 12, base_family = "") 
{
    theme_noGrid(base_size = base_size, base_family = base_family) %+replace% 
        theme(
                panel.grid.major.y = element_line(colour = "grey78",size = 0.3),
                panel.grid.minor.y = element_line(colour = "grey88",size = 0.25),
            )
}

theme_grid_x <- function (base_size = 12, base_family = "") 
{
    theme_noGrid(base_size = base_size, base_family = base_family) %+replace% 
        theme(
                panel.grid.major.x = element_line(colour = "grey78",size = 0.3),
                panel.grid.minor.x = element_line(colour = "grey88",size = 0.25),
            )
}


theme_fullBlank <- function (base_size = 12, base_family = "") 
{
    theme_grey(base_size = base_size, base_family = base_family) %+replace% 
        theme(
                axis.line = element_blank(),        
                axis.text = element_blank(), 
                axis.ticks = element_blank(), 
                axis.title = element_blank(),               
                legend.key = element_rect(colour = "grey80"), 
                panel.background = element_rect(fill = "white", colour = NA), 
                panel.border = element_rect(fill = NA,colour = NA),
                panel.grid.major = element_blank(),                  
                panel.grid.minor = element_blank(), 
                strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2)
            )
}

# Functions to create a vertical dodge function
# Detect and prevent collisions.
# Powers dodging, stacking and filling.
collidev <- function(data, height = NULL, name, strategy, check.height = TRUE) {
    # Determine height
    if (!is.null(height)) {
        # height set manually
        if (!(all(c("ymin", "ymax") %in% names(data)))) {
            data$ymin <- data$y - height / 2
            data$ymax <- data$y + height / 2
        }
    } else {
        if (!(all(c("ymin", "ymax") %in% names(data)))) {
            data$ymin <- data$y
            data$ymax <- data$y
        }
        
        # height determined from data, must be floating point constant
        heights <- unique(data$ymax - data$ymin)
        heights <- heights[!is.na(heights)]
        
        #   # Suppress warning message since it's not reliable
        #     if (!zero_range(range(heights))) {
        #       warning(name, " requires constant height: output may be incorrect",
        #         call. = FALSE)
        #     }
        height <- heights[1]
    }
    
    # Reorder by x position, relying on stable sort to preserve existing
    # ordering, which may be by group or order.
    data <- data[order(data$ymin), ]
    
    # Check for overlap
    intervals <- as.numeric(t(unique(data[c("ymin", "ymax")])))
    intervals <- intervals[!is.na(intervals)]
    
    if (length(unique(intervals)) > 1 & any(diff(scale(intervals)) < -1e-6)) {
        warning(name, " requires non-overlapping y intervals", call. = FALSE)
        # This is where the algorithm from [L. Wilkinson. Dot plots.
        # The American Statistician, 1999.] should be used
    }
    
    if (!is.null(data$xmax)) {
        plyr::ddply(data, "ymin", strategy, height = height)
    } else if (!is.null(data$x)) {
        data$xmax <- data$x
        data <- plyr::ddply(data, "ymin", strategy, height = height)
        data$x <- data$xmax
        data
    } else {
        stop("Neither x nor xmax defined")
    }
}

# Stack overlapping intervals.
# Assumes that each set has the same horizontal position
pos_stackv <- function(df, height) {
    if (nrow(df) == 1) return(df)
    
    n <- nrow(df) + 1
    x <- ifelse(is.na(df$x), 0, df$x)
    if (all(is.na(df$y))) {
        heights <- rep(NA, n)
    } else {
        heights <- c(0, cumsum(x))
    }
    
    df$xmin <- heights[-n]
    df$xmax <- heights[-1]
    df$x <- df$xmax
    df
}

# Stack overlapping intervals and set height to 1.
# Assumes that each set has the same horizontal position.
pos_fillv <- function(df, height) {
    stacked <- pos_stackv(df, height)
    stacked$xmin <- stacked$xmin / max(stacked$xmax)
    stacked$xmax <- stacked$xmax / max(stacked$xmax)
    stacked$x <- stacked$xmax
    stacked
}

# Dodge overlapping interval.
# Assumes that each set has the same horizontal position.
pos_dodgev <- function(df, height) {
    n <- length(unique(df$group))
    if (n == 1) return(df)
    
    if (!all(c("ymin", "ymax") %in% names(df))) {
        df$ymin <- df$y
        df$ymax <- df$y
    }
    
    d_height <- max(df$ymax - df$ymin)
    
    # df <- data.frame(n = c(2:5, 10, 26), div = c(4, 3, 2.666666,  2.5, 2.2, 2.1))
    # ggplot(df, aes(n, div)) + geom_point()
    
    # Have a new group index from 1 to number of groups.
    # This might be needed if the group numbers in this set don't include all of 1:n
    groupidy <- match(df$group, sort(unique(df$group)))
    
    # Find the center for each group, then use that to calculate xmin and xmax
    df$y <- df$y + height * ((groupidy - 0.5) / n - .5)
    df$ymin <- df$y - d_height / n / 2
    df$ymax <- df$y + d_height / n / 2
    
    df
}


#' Adjust position by dodging overlaps to the side.
#'
#' @inheritParams ggplot2::position_identity
#' @param height Dodging height, when different to the height of the individual
#'   elements. This is useful when you want to align narrow geoms with wider
#'   geoms. See the examples for a use case.
#' @family position adjustments
#' @export
position_dodgev <- function(height = NULL) {
    ggproto(NULL, PositionDodgeV, height = height)
}



PositionDodgeV <- ggproto("PositionDodgeV", Position,
                          required_aes = "y",
                          height = NULL,
                          setup_params = function(self, data) {
                              if (is.null(data$ymin) && is.null(data$ymax) && is.null(self$height)) {
                                  warning("height not defined. Set with `position_dodgev(height = ?)`",
                                          call. = FALSE)
                              }
                              list(height = self$height)
                          },
                          
                          compute_panel = function(data, params, scales) {
                              collidev(data, params$height, "position_dodgev", pos_dodgev, check.height = FALSE)
                          }
)

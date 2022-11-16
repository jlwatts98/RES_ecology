##### Functions #####

firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

ggbiplot2 <- function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                       obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                       ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                       alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                       varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                       color = muted("red"), # <- add new arguments to the function
                       linetype = "solid",
                       alpha_arrow = 1)
  
g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, xend = xvar, yend = yvar),
                        arrow = arrow(length = unit(1/2, "picas")),
                        color = color, linetype = linetype, alpha = alpha_arrow)


##############################################
#
# General description:
#
#   The following function plots a set of genomic ranges
#   

# Input:
#
#     Ranges:     ranges to plot
#     ColV:       vector of colors (as long as GRanges)


# Output:
#   

##############################################

PlotGRanges <- function(Ranges, ColV = rep("grey", length(Ranges))){
  
  if (length(Ranges) != length(ColV)){
    warning("Ranges and color vector don't have the same length!\n")
  }
  plot(c(min(start(Ranges)), max(end(Ranges))), c(1, 1.1 * length(Ranges)), 
       yaxt = "n", xlab = "Genomic coordinate", ylab = "", type = "n",
       bty = "n")
  rect(start(Ranges), seq_along(Ranges), end(Ranges), seq_along(Ranges) + 0.3,
       col = ColV)
 }



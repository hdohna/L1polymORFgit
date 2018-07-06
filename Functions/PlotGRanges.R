##############################################
#
# General description:
#
#   The following function plots a set of genomic ranges
#   

# Input:
#
#     Ranges:     ranges to plot


# Output:
#   

##############################################

PlotGRanges <- function(Ranges){
  
  plot(c(min(start(Ranges)), max(end(Ranges))), c(1, length(Ranges)), 
       yaxt = "n", xlab = "Genomic coordinate", ylab = "", type = "n",
       bty = "n")
  rect(start(Ranges), seq_along(Ranges), end(Ranges), seq_along(Ranges) + 0.3,
       col = "grey")
 }



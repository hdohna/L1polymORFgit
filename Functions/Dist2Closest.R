##############################################
#
# General description:
#
#   The following function calculates the distances from one set of genomic
#   ranges to the closest range from a second set of genomic ranges
#   

# Input:
#
#     GR1: first set of genomic ranges
#     GR2: first set of genomic ranges

# Output:
#   
#     ...

##############################################

Dist2Closest <- function(GR1, GR2){
  DistObj <- distanceToNearest(GR1, GR2, ignore.strand = T) 
  Dists <- DistObj@elementMetadata@listData$distance
  
  # Return a warning message if not all GR1s received a distance
  if (length(Dists) != length(GR1)) {
    warning("Not the same number of distances as input ranges!\n")
  }
  Dists
}
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
  DistObj <- distanceToNearest(GR1, GR2, ignore.strand = T, select="all") 
  Dists <- DistObj@elementMetadata@listData$distance
  
  # Return a warning message if not all GR1s received a distance
  if (length(Dists) != length(GR1)) {
    Dists <- rep(NA, length(GR1))
    Dists[DistObj@from] <- DistObj@elementMetadata@listData$distance
    warning(sum(is.na(Dists)), " elements of GR1 do not have a distance!\n")
  }
  Dists
}
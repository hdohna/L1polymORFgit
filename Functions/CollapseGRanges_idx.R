##############################################
#
# General description:
#
#   The following function assigns IDs to ranges and gives overlapping
#   ranges the same ID
#   

# Input:
#
#     GR: ranges to collapse
#     NewSize: integer value for new range size
#     blnPairwise: whether collapse by pairwise overlap


# Output:
#   
#    ID for each range (overlapping ranges get same ID)

# Comments:
#   
#    Requires package GenomicRanges. This function is used in CollapseClosePos_idx

##############################################

CollapseGRanges_idx <- function(GR, NewSize = NULL, blnPairwise = F){
  
  if(! is.null(NewSize)){
    GR <- resize(GR, NewSize, fix = "center")
  }
  
  if (blnPairwise) {
    OL <- findOverlaps(GR, GR)
    
    # Assign to all overlapping ranges the smallest index
    idx <- sapply(1:length(GR), function(x){
      blnOL <- OL@from ==x | OL@to == x
      min(c(OL@from[blnOL], OL@to[blnOL]))
    })
    
  } else {
    GRUnion <- union(GR, GR)
    OL  <- findOverlaps(GR, GRUnion)
    if(length(OL@from) != length(GR)) {
      warning("No 1:1 match between original and collapsed ranges!\n")
    }
    idx <- OL@to
  }
  
  
}



##############################################
#
# General description:
#
#   The following function creates a unique set of genomic ranges by
#   retaining for each set of overlapping ranges the longest one
#   

# Input:
#
#     Ranges: ranges to analyze
#     Group: grouping variable

# Output: 
#    NewRanges: New ranges with overlaps removed
#    

##############################################

UniqueGRanges <- function(Ranges, Group = NULL){
  
  # If no grouping variable is supplied
  if (is.null(Group)){
    RangesOL  <- findOverlaps(Ranges, Ranges)
    idxDupl   <- unique(RangesOL@from[duplicated(RangesOL@from)])
    idxAppend <- sapply(idxDupl, function(x) {
      idxSet <- RangesOL@from[RangesOL@from == x]
      idxSet[which.max(width(Ranges[idxSet]))]
    })
    NewRanges <- c(Ranges[-idxDupl], Ranges[idxAppend]) 
  
  # If grouping variable is supplied
  } else {
    if (length(Ranges) != length(Group)){
      stop("Ranges and Group must be same length")
    }
    GroupUnique <- unique(Group)
    idxRetain <- sapply(GroupUnique, function(x) {
      idxSet <- which(Group == x)
      idxSet[which.max(width(Ranges[idxSet]))]
    })
    NewRanges <- Ranges[idxRetain]
  }
}



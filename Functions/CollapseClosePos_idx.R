##############################################
#
# General description:
#
#   The following function assigns IDs to genomic positions and assigns to
#   close positions the same ID
#   

# Input:
#
#     DF:          input data frame
#     ChromCol:    name of column containing chromosome
#     PosCol:      name of column containing position
#     OLRange:     integer value for range so that positions within that range
#                  get the same ID
#     blnPairwise: whether collapse by pairwise overlap


# Output:
#   
#    ID for each position (close positions get same ID)

# Comments:
#   
#    Requires package GenomicRanges and function CollapseGRanges_idx

##############################################

CollapseClosePos_idx <- function(DF, ChromCol, PosCol, OLRange = 10, 
                                 blnPairwise = F){
  
  # Create genomic ranges from positions
  GR <- makeGRangesFromDataFrame(DF, seqnames.field   = ChromCol,
                                          start.field = PosCol,
                                          end.field   = PosCol)
  
  # Get indices for GRs based on overlaps
  CollapseGRanges_idx(GR, NewSize = 2*OLRange, blnPairwise = blnPairwise)
}



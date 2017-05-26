##############################################
#
# General description:
#
#   The following function takes a read list (as created by the function 
#   scanBam) and identifies regions that are covered by the minimum number
#   of ZMW IDs (PacBio) and minimum number of reads

# Input:
#
#     RL: list with read info (as created by the function scanBam)
#     MinZMWCover: minimum number of different ZMW ID required (usually 2)
#     MinReadCover: minimum number of reads.
#     ChrLen: Chromosome length (can be any large number extending beyond 
#        covered region)

# Output:
#   
#    MinCoverGR: Genomic ranges with minimum coverage

# Comments:
#   
#    Requires function ZMWCoverage and ReadList2GRanges

##############################################

MinCoverRangesPacBio <- function(RL, MinZMWCover = 2, 
                                 MinReadCover = 3, ChrLen = 10^6){
  
  # Get coverage of unique ZMWs
  ZCover <- ZMWCoverage(RL, ChrLen)
  
  # Get read coverage
  RGR    <- ReadList2GRanges(RL)
  RCover <- coverage(RGR)[RL$rname[1]][[1]]
  
  # Get regions with minimum coverage
  ZCoverIR <- slice(ZCover, lower = MinZMWCover, rangesOnly = T)
  RCoverIR <- slice(RCover, lower = MinReadCover, rangesOnly = T)
  ZCoverGR <- GRanges(RL$rname[1], ZCoverIR)
  RCoverGR <- GRanges(RL$rname[1], RCoverIR)
  
  # Intersect regions with minimum coverage
  MinCoverGR <- intersect(ZCoverGR, RCoverGR)
}



##############################################
#
# General description:
#
#   The following function takes a read list (as created by the function 
#   scanBam), a genomic range and gets the consensus sequence per ZMW

# Input:
#
#     RL: list with read info (as created by the function scanBam)
#     GR: Genomic range within which to call SNPs
#     IncludeStar: option to include stars. Reads that do not extend over the
#        entire range of GR are extended with stars. Consensus sequences can
#        include the stars.

# Output:
#   
#    SNPpos: SNP positions

# Comments:
#   
#    Requires function SeqMatFromReads, SeqFromCigar and ReadList2GRanges. 
#    RL cannot contain NA. Sequences are extened using *

##############################################

ConsensusFromReads <- function(RL, GR, IncludeStar = F){
  
  # Get a matrix of sequences
  SeqMat <- SeqMatFromReads(RL, GR)
  
  # Get the consensus sequence
  apply(SeqMat, 1, FUN = function(x) {
    NucCount <- table(x)
    NucCount <- NucCount[names(NucCount) != "*" | IncludeStar]
    names(NucCount)[which.max(NucCount)]
  })
}



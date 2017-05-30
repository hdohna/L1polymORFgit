##############################################
#
# General description:
#
#   The following function takes a read list (as created by the function 
#   scanBam), a genomic range and gets a matrix of phred scores

# Input:
#
#     RL: list with read info (as created by the function scanBam)
#     GR: Genomic range within which to call SNPs

# Output:
#   
#    Matrix of numeric phred values(row = sequence position, column = sequence)

# Comments:
#   
#    Requires function SeqFromCigar and ReadList2GRanges. 
#    RL cannot contain NA. Sequences are extened using *

##############################################

PhredMatFromReads <- function(RL, GR){
  
  # Turn reads into genomic ranges and only retain the ones overlapping with GR
  RGR         <- ReadList2GRanges(RL)
  StartAll    <- start(GR)
  EndAll      <- end(GR)
  blnInRange  <- overlapsAny(RGR, GR)
  RL          <- lapply(RL, function(x) x[blnInRange])
  
  # Get a matrix of phred scores
  sapply(seq_along(RL$pos), function(j) {
    PhredV      <- PhredFromCigar(RL$cigar[j], RL$qual[j])
    PhredStart  <- max(1, StartAll - RL$pos[j] + 1)
    PhredEnd    <- min(length(PhredV), EndAll - RL$pos[j] + 1)
    NrPrepend <- max(0, RL$pos[j] - StartAll)
    Prepend   <- rep(-1, NrPrepend)
    NrAppend  <- max(0, EndAll - RL$pos[j] - PhredEnd + 1)
    Append    <- rep(-1, NrAppend)
    c(Prepend, PhredV[PhredStart:PhredEnd], Append)
  })
}



##############################################
#
# General description:
#
#   The following function takes a read list (as created by the function 
#   scanBam), a genomic range and gets a matrix of of sequences

# Input:
#
#     RL: list with read info (as created by the function scanBam)
#     GR: Genomic range within which to call SNPs
#     MinProp: minimum proportion of reads that should contain SNP

# Output:
#   
#    Matrix of sequences (row = sequence position, column = sequence)

# Comments:
#   
#    Requires function SeqFromCigar and ReadList2GRanges. 
#    RL cannot contain NA. Sequences are extened using *

##############################################

SeqMatFromReads <- function(RL, GR){
  
  # Turn reads into genomic ranges and only retain the ones overlapping with GR
  RGR         <- ReadList2GRanges(RL)
  StartAll    <- start(GR)
  EndAll      <- end(GR)
  blnInRange  <- overlapsAny(RGR, GR)
  RL          <- lapply(RL, function(x) x[blnInRange])
  
  # Get a matrix of sequences
  sapply(seq_along(RL$pos), function(j) {
    SeqV      <- SeqFromCigar(RL$cigar[j], RL$seq[j])
    SeqStart  <- max(1, StartAll - RL$pos[j] + 1)
    SeqEnd    <- min(length(SeqV), EndAll - RL$pos[j] + 1)
    NrPrepend <- max(0, RL$pos[j] - StartAll)
    Prepend   <- rep("*", NrPrepend)
    NrAppend  <- max(0, EndAll - RL$pos[j] - SeqEnd + 1)
    Append    <- rep("*", NrAppend)
    c(Prepend, SeqV[SeqStart:SeqEnd], Append)
  })
}



##############################################
#
# General description:
#
#   The following function converts a read list (as created by the function 
#   scanBam) into genomic ranges (1 per read).

# Input:
#
#     RL: list with read info (as created by the function scanBam)

# Output:
#   
#    : ...

# Comments:
#   
#    Requires packages csaw

##############################################

ReadList2GRanges <- function(RL){
  
  # Get start positions, end positions 
  Starts <- RL$pos
  Ends   <- Starts + sapply(RL$cigar, ReadLengthFromCigar, USE.NAMES = F)
  
  # Create genomic ranges
  GRanges(seqnames = RL$rname, ranges = IRanges(start = Starts, end = Ends))
}



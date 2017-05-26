##############################################
#
# General description:
#
#   The following function reads a read list (as created by the function 
#   scanBam) and creates an Rle object with the coverage of different
#   ZMW IDs.

# Input:
#
#     RL: list with read info (as created by the function scanBam)
#     ChrLength: length of chromosome

# Output:
#   
#    : ...

# Comments:
#   
#    Requires packages csaw

##############################################

ZMWCoverage <- function(RL, ChrLength){
  
  # Get start positions, end positions and ZMW IDs
  ZMW_ID <- sapply(RL$qname, function(x) strsplit(x, "/")[[1]][2], USE.NAMES = F)
  Starts <- RL$pos
  Ends   <- Starts + sapply(RL$cigar, ReadLengthFromCigar, USE.NAMES = F)
  
  # Get three ordered vectors
  ZMWV       <- c(ZMW_ID, ZMW_ID)
  PosV       <- c(Starts, Ends)
  ChangeV    <- c(rep(1, length(Starts)), rep(-1, length(Ends)))
  POrder     <- order(PosV)
  PosV       <- PosV[POrder]
  ChangeV    <- ChangeV[POrder]
  ZMWV       <- ZMWV[POrder]
  ZMWVNoDupl <- ZMWV[!duplicated(ZMWV)]
  
  # Initialize loop to create coverage RLE
  RleValues       <- 0
  RleLengths      <- NULL
  ChangePoint     <- 1
  ZMWCount        <- rep(0, length(ZMWVNoDupl))
  names(ZMWCount) <- ZMWVNoDupl
  CurrentCover    <- 0
  
  # Loop over points and increment RLE values
  for (i in 1:length(PosV)){
    ZMWCount[ZMWV[i]] <- ZMWCount[ZMWV[i]] + ChangeV[i]
    NewCover          <- sum(ZMWCount > 0)
    if (NewCover != CurrentCover){
      RleValues   <- c(RleValues, NewCover)
      RleLengths  <- c(RleLengths, PosV[i] - ChangePoint)
      ChangePoint <- PosV[i]
      CurrentCover <- NewCover
    }
  }
  
  # Add last run length value
  RleLengths <- c(RleLengths, ChrLength - ChangePoint)
  Rle(RleValues, RleLengths)
}



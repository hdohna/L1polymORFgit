##############################################
#
# General description:
#
#   The following function reads in a fastq file stepwise and collects
#   a vector or read lengths
#   

# Input:
#
#     FilePath: (character) path to fastq file that will read
#     NrLines2Read: (integer) number of lines to read per iteration

# Output:
# 
#    ReadLengths: (integer) vector of read lengths 
#     ...

##############################################

getReadLengthFromFastQ <- function(FilePath,
                         NrLines2Read = 10^5){
  
  cat("***********************************************************\n")
  cat("**                                                       **\n")
  cat("**    Running function getReadLengthFromFastQ ...          **\n")
  cat("**                                                       **\n")
  cat("***********************************************************\n")
  
  ReadLengths  <- c()
  LinesSkip    <- 0
  NrLinesRead  <- NrLines2Read
  while(NrLinesRead >= NrLines2Read){
    Lines <- scan(FilePath, skip = LinesSkip, nlines = NrLines2Read,
                  what = "character")
    idxAt       <- which(substr(Lines, 1, 1) == "@")
    ReadLengths <- c(ReadLengths, nchar(Lines[idxAt + 1]))
    NrLinesRead <- length(Lines)
    LinesSkip   <- LinesSkip + NrLinesRead
    cat("Processed", LinesSkip, "lines\n")
  }
  ReadLengths
}
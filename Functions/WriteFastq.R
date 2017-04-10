# General description:
#
#    This function writes a vector of reads out as a FastQ file and writes a 
#    sample file to be used by qAlign

# Arguments:
#   
#    Reads: character vector of read sequences
#    ReadNames: character vector of read names
#    FilePath: character string indicating file path of fastq file
#        to be written out

# Output:
#   
#    fastq file saved at location specified in FilenamePath

# Comment:

WriteFastq <- function(Reads, 
                       ReadNames = NULL, ReadQual = NULL,
                       FilePath){
  
  # Create read names if not provided
  if (is.null(ReadNames)){
    ReadNames <- paste("@", 1:length(Reads), sep = "")
  }
  
  # Add @ to read names if not present
  RReadNames1stChar <- substr(ReadNames, 1, 1)
  blnNotAt <- RReadNames1stChar != "@"
  ReadNames[blnNotAt] <- paste("@", ReadNames[blnNotAt], sep = "")
  
  # Check that vector of reads and read names are the same length
  if (length(Reads) != length(ReadNames)){
    stop("Reads and read names must be the same length!")
  }
  
  # Create read quality if not provided
  if (is.null(ReadQual)){
    ReadQual <- sapply(nchar(Reads), function(x) {
      QVect <- rep("~", x)
      paste(QVect, collapse = "")
    })
  }
  
  
  # Create lines of the fastq file and save them
  FastQLines3 <- rep("+", length(Reads))
  FastQLines4 <- ReadQual
  FastQLinesAll <- rep(NA, 4 * length(Reads))
  FastQLinesAll[seq(1, length(FastQLinesAll) - 3, 4)] <- ReadNames
  FastQLinesAll[seq(2, length(FastQLinesAll) - 2, 4)] <- Reads
  FastQLinesAll[seq(3, length(FastQLinesAll) - 1, 4)] <- FastQLines3
  FastQLinesAll[seq(4, length(FastQLinesAll) - 0, 4)] <- FastQLines4
  writeLines(FastQLinesAll, FilePath)
  
}



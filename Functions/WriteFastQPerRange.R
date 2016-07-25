##############################################
#
# General description:
#
#     The following function takes a set of genomic ranges, reads all read IDs
#     in that range from a bam file, get's matching IDs from a fastq file and
#     then saves each set of reads to a separate fastq file. This is meant for
#     paired-end reads

# Input:
#
#     Ranges: Genomic ranges. A separate fastq file will be saved for each 
#          range
#     InBamfilePath: character string indicating path to bam file 
#     InFastQfilePaths: character vector of length 2 indicating paths to the 
#          two fastq files containing the paired reads
#     OutFilePaths: vector of character strings indicating paths to all the
#          fastq files that are written out. Should be same length as Ranges
#     NrReadsPerIter: integer value indicating how many reads from the fastq
#          file should be read in per iteration
#     DefaultWriteMode: single character equal to either 'w' or 'a' to write
#          to a new file or append to an existing file, respectively.


# Output:
#   
#    : ...

# Comments:
#    This function requires packages Rsamtools and ShortRead

##############################################


WriteFastQPerRange <- function(Ranges, InBamfilePath, InFastQfilePaths,
  OutFilePaths, NrReadsPerIter = 10^6, DefaultWriteMode = 'w') {
  
  cat("****************************************************\n")
  cat("**                                                **\n")
  cat("**    Running function WriteFastQPerRange ...     **\n")
  cat("**                                                **\n")
  cat("****************************************************\n\n")

    # Check that Ranges and OutFilePaths have the same length
  if (length(Ranges) !=  length(OutFilePaths)) {
    stop("Ranges and OutFilePaths have to have the same length! \n")
  }
  
  # Get read IDs per range
  param <- ScanBamParam(which = Ranges, what = "qname")
  ReadIDsPerRange <- scanBam(file = InBamfilePath, param = param)
  
  # Create a vector of all read names and range indices
  AllReadNames    <- unlist(ReadIDsPerRange) 
  NrReadsPerRange <- sapply(ReadIDsPerRange, function(x) length(x$qname))
  AllReadIndices  <- rep(1:length(ReadIDsPerRange), NrReadsPerRange) 

  # Open one connection per fastq file and intitialize counter variables
  Stream1 <- FastqStreamer(InFastQfilePaths[1], NrReadsPerIter)
  if (length(InFastQfilePaths) > 1){
    Stream2 <- FastqStreamer(InFastQfilePaths[2], NrReadsPerIter)
  }
  NrReadsRead  <- 1
  ReadCounter  <- 0
  WriteCounter <- 0
  wqModes      <- rep(DefaultWriteMode, length(OutFilePaths))
  
  # Loop through fastq file and append to little range-specific fastq files
  while (NrReadsRead > 0){
    Reads1  <- yield(Stream1)
    if (length(InFastQfilePaths) > 1){
      Reads2  <- yield(Stream2)
    }
    NrReadsRead  <- length(Reads1)
    ReadCounter  <- ReadCounter + NrReadsRead
    ReadIDChar   <- as.character(Reads1@id)
    ReadIDsuffix <- substr(ReadIDChar, nchar(ReadIDChar) - 7, nchar(ReadIDChar))
    ReadIDs      <- substr(ReadIDChar, 1, nchar(ReadIDChar) - 8)
    ReadSubset   <- ReadIDs %in% AllReadNames
    ReadIDs      <- ReadIDs[ReadSubset]
    Reads1       <- Reads1[ReadSubset]
    if (length(InFastQfilePaths) > 1){
      Reads2       <- Reads2[ReadSubset]
    }
    cat("Total of", ReadCounter, "reads read \n")
    # Loop over ranges and create a fastq file per range
    cat("Looping over ranges and writing to fastq files per range\n")
    browser()
    if (sum(ReadSubset) > 0) {
      Indices <- unique(AllReadIndices[AllReadNames %in% ReadIDs])
      for (i in Indices){
        
        # Get names of reads in current range, subset reads of current chunk and 
        # append to fastq files
        ReadNames   <- ReadIDsPerRange[[i]]$qname
        ReadSubset  <- ReadIDs %in% ReadNames
        Reads1Local <- Reads1[ReadSubset]
         writeFastq(Reads1Local, OutFilePaths[i], mode = wqModes[i], 
                   compress = F, full = T) 
        wqModes[i] <- "a"
        if (length(InFastQfilePaths) > 1){
          Reads2Local <- Reads2[ReadSubset]
          writeFastq(Reads2Local, OutFilePaths[i], mode = wqModes[i], 
                   compress = F, full = T) 
        }
        WriteCounter <- WriteCounter + 2*sum(ReadSubset)
        
      }
    }
    cat("Total of", WriteCounter, "written out \n")
  }
  close(Stream1)
  if (length(InFastQfilePaths) > 1){
    close(Stream2)
  }
  
}




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
#     OutFilePaths: vector of character strings indicating paths to all the
#          fastq files that are written out. Should be same length as Ranges
#     NrReadsPerIter: integer value indicating how many reads from the fastq
#          file should be read in per iteration
#     DefaultWriteMode: single character equal to either 'w' or 'a' to write
#          to a new file or append to an existing file, respectively.
#     IdChar2Remove: number of characters to remove from fastq file read IDs
#          so that they mtch the bam file read IDs


# Output:
#   
#    : ...

# Comments:
#    This function requires packages Rsamtools and ShortRead

##############################################


WriteFastQPerRangeFromBam <- function(Ranges, InBamfilePath, 
  OutFilePaths) {
  
  cat("***********************************************************\n")
  cat("**                                                       **\n")
  cat("**    Running function WriteFastQPerRangeFromBam ...     **\n")
  cat("**                                                       **\n")
  cat("***********************************************************\n")
  
  
  # Check that Ranges and OutFilePaths have the same length
  if (length(Ranges) !=  length(OutFilePaths)) {
    stop("Ranges and OutFilePaths have to have the same length! \n")
  }
  
  # Status messages
  cat("Writing out fastq files from bam file\n")
  cat("Number of fastq files to write out:", length(Ranges), "\n")
    
  # Get info per range
  param <- ScanBamParam(which = Ranges, what = scanBamWhat())
  ReadsPerRange <- scanBam(file = InBamfilePath, param = param)
    
  # Loop over ranges and write out fastq file per range
  for(i in 1:length(ReadsPerRange)){
      RL <- ReadsPerRange[[i]]
      blnNegStrand <- RL$strand == "-"
      RL$seq[blnNegStrand] <- reverseComplement(RL$seq[blnNegStrand])
      RL$qual[blnNegStrand] <- reverse(RL$qual[blnNegStrand])
      cat("Writing out file", i, "out of", length(ReadsPerRange), "\n")
      WriteFastq(Reads = as.character(RL$seq), 
                 ReadNames = RL$qname, 
                 ReadQual = as.character(RL$qual),
                 FilePath = OutFilePaths[i])
      cat("File written out:", OutFilePaths[i], "\n")
  }
}




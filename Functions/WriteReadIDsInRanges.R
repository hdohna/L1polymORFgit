##############################################
#
# General description:
#
#     The following function takes a set of genomic ranges, reads all read IDs
#     in thes ranges from a bam file, and writes all the IDs in those ranges
#     to a text file

# Input:
#
#     Ranges: Genomic ranges. A separate fastq file will be saved for each 
#          range
#     InBamfilePath: character string indicating path to bam file 
#     OutFilePath: character string indicating paths to the text file that 
#          contains all the read IDs


# Output:
#   
#    : ...

# Comments:
#    This function requires package Rsamtools

##############################################


WriteReadIDsInRanges <- function(Ranges, InBamfilePath, 
  OutFilePath) {
  
  # Get read IDs per range
  param <- ScanBamParam(which = Ranges, what = "qname")
  ReadIDsPerRange <- scanBam(file = InBamfilePath, param = param)
  
  # Create a vector of all unique read names
  ReadIDsPerRange    <- unlist(ReadIDsPerRange) 
  ReadIDsPerRange    <- unique(ReadIDsPerRange) 
  
  # Write all read IDs to a text file 
  writeLines(ReadIDsPerRange, OutFilePath)
}




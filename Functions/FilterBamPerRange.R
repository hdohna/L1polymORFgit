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
#     OutBamFilePaths: vector of character strings indicating paths to all the
#          bam files that are written out. Should be same length as Ranges
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


FilterBamPerRange <- function(Ranges, InBamfilePath,
  OutBamFilePaths) {
  
  cat("****************************************************\n")
  cat("**                                                **\n")
  cat("**    Running function FilterBamPerRangeByID ...     **\n")
  cat("**                                                **\n")
  cat("****************************************************\n\n")

  # Check that Ranges and OutBamFilePaths have the same length
  if (length(Ranges) !=  length(OutBamFilePaths)) {
    stop("Ranges and OutBamFilePaths have to have the same length! \n")
  }
  
  # Filter per range
    for (i in 1:length(Ranges)){
    cat("Filtering reads of range", i, "of", length(Ranges), "\n")
    filterBam(InBamfilePath, OutBamFilePaths[i], 
              param=ScanBamParam(what=scanBamWhat(), which = Ranges[i]))
    }
  OutBamFilePaths
}




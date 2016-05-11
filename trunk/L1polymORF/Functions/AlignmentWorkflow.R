##############################################
#
# General description:
#
#   The following function starts with a fastq file name, constructs and calls
#   commands in the workflow from alignment, de-duplication, idexing and sorting

# Input:
#
#  Folders and file names:
#     FastqFile:     name of fastq file
#     ReferenceFile: name of reference file
#     SamFile:       name of sam file produced by aligning fastq to reference
#     BamFileSorted:   name of bam file after sorting
#     BamFileDedup:  name of bam file after deduplication

#  Run parameters:
#     blnAlign: boolean indicator whether alignment step should be performed
#     blnDedup: boolean indicator whether deduplication step should be performed
#     blnSort: boolean indicator whether sorting step should be performed

#  Run commands:
#     AlignCommand: character string with command to run alignment 
#     AddGroupCmd: character string with command to add read groups to bam 
#         files
#     AddGroupOptions: character vector specifying options add read groups to
#         bam files
#     DedupCommand: character string with command to run depuplication 
#     SortCommand: character string with command to run sorting 

# Output:
#   
#    Bam file after preprocessing

##############################################

AlignmentWorkflow <- function(FastqFile, 
                              ReferenceFile, 
                              SamFile = NULL,
                              BamFile = NULL,
                              BamFileSorted = NULL,
                              BamFileDedup = NULL,
                              MetricsFileDedup = NULL,
                              BamFileUnique = NULL,
                              BamFileSorted2 = NULL,
                              blnAlign = F,
                              blnSam2SortedBam = F,
                              blnDedup = F,
                              blnRemMapQ0 = F,
                              blnSort = F,
                              AlignCommand = '/home/txw/bwa/bwa-0.7.12/bwa mem',
                              DedupCommand ='java -jar /home/txw/picard/picard-tools-1.131/picard.jar MarkDuplicates',
                              DedupOptions = c('REMOVE_DUPLICATES=true ASSUME_SORTED=true', 
                               'TMP_DIR=/home/hzudohna/tmp', 
                               'VALIDATION_STRINGENCY=SILENT', 
                               'MAX_RECORDS_IN_RAM=2000000', 
                               'PROGRAM_RECORD_ID=MarkDuplicates', 
                               'PROGRAM_GROUP_NAME=MarkDuplicates', 
                               'MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000', 
                               'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000', 
                               'SORTING_COLLECTION_SIZE_RATIO=0.25', 
                               'READ_NAME_REGEX=[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*', 
                               'OPTICAL_DUPLICATE_PIXEL_DISTANCE=100', 
                               'VERBOSITY=INFO', 
                               'QUIET=false', 
                               'COMPRESSION_LEVEL=5', 
                               'CREATE_INDEX=false', 
                               'CREATE_MD5_FILE=false')
){
  
  # Create file names
  if (is.null(SamFile)) {
    FilePrefix <- strsplit(FastqFile, ".fastq")[[1]][1]
    SamFile    <- paste(FilePrefix, ".sam", sep = "")
  }
  if (is.null(BamFile)) {
    FilePrefix <- strsplit(SamFile, ".sam")[[1]][1]
    BamFile    <- paste(FilePrefix, ".bam", sep = "")
  }
  if (is.null(BamFileSorted)) {
    FilePrefix  <- strsplit(BamFile, ".bam")[[1]][1]
    BamFileSorted <- paste(FilePrefix, ".sorted.bam", sep = "")
  }
  if (is.null(MetricsFileDedup)) {
    FilePrefix   <- strsplit(BamFileSorted, ".sorted.bam")[[1]][1]
    MetricsFileDedup <- paste(FilePrefix, ".dedup.metrics", sep = "")
  }
  if (is.null(BamFileDedup)) {
    FilePrefix   <- strsplit(BamFileSorted, ".sorted.bam")[[1]][1]
    BamFileDedup <- paste(FilePrefix, ".dedup.bam", sep = "")
  }
  if (is.null(BamFileUnique)) {
    FilePrefix    <- strsplit(BamFileDedup, ".bam")[[1]][1]
    BamFileUnique <- paste(FilePrefix, ".unique.bam", sep = "")
  }
  if (is.null(BamFileSorted2)) {
    FilePrefix  <- strsplit(BamFileUnique, ".bam")[[1]][1]
    BamFileSorted2 <- paste(FilePrefix, ".sorted.bam", sep = "")
  }
  
  
  # Run alignment
  if (blnAlign){
    cat("\n\n*******************************************************\n")
    cat("**                                                   **\n")
    cat("**    Aligning", FastqFile," ...             **\n")
    cat("**                                                   **\n")
    cat("*******************************************************\n\n")
    
    CmdLine <- paste(AlignCommand,  ReferenceFile, FastqFile, ">", SamFile)
    system(CmdLine)
  }
  
  # Turn sam file into sorted and indexed bam file
  if (blnSam2SortedBam){
    cat("\n\n*******************************************************\n")
    cat("**                                                   **\n")
    cat("**    Converting", SamFile," to bam file   **\n")
    cat("**                                                   **\n")
    cat("*******************************************************\n\n")
    
    # Command line file for conversion
    CmdLine <- paste("/home/txw/samtools/samtools-1.2/samtools view", 
                     SamFile, "-b -h -o", BamFile)
    system(CmdLine)
    
    cat("\n\n*******************************************************\n")
    cat("**                                                   **\n")
    cat("**    Sorting", BamFile,"    **\n")
    cat("**                                                   **\n")
    cat("*******************************************************\n\n")
    
    # Command line file for indexing
    CmdLine <- paste("/home/txw/samtools/samtools-1.2/samtools index", 
                     BamFile)
    system(CmdLine)
    
    # Command line file for sorting
    CmdLine <- paste("/home/txw/samtools/samtools-1.2/samtools sort -o", 
                     BamFileSorted, "-T", paste(BamFile, ".tmp", sep = ""), 
                     BamFile)
    system(CmdLine)
    
  }

  # Run deduplication
  if (blnDedup){
    cat("\n\n*******************************************************\n")
    cat("**                                                   **\n")
    cat("**    Deduplicating", BamFileSorted," ...             **\n")
    cat("**                                                   **\n")
    cat("*******************************************************\n\n")
    
    # Create command lines
    DedupInputLine   <- paste('INPUT=', BamFileSorted, sep = "")
    DedupOutputLine  <- paste('OUTPUT=', BamFileDedup, sep = "")
    DedupMetricsLine <- paste('METRICS_FILE=', MetricsFileDedup, sep = "")

    CmdLine <- paste(DedupCommand, DedupInputLine, DedupOutputLine, 
                     DedupMetricsLine, DedupOptions)
    system(CmdLine)
  }
  
  # Run removal of zero mapping quality
  if (blnRemMapQ0){
    cat("\n\n*******************************************************\n")
    cat("**                                                   **\n")
    cat("**    Remove mapq 0 from", BamFileDedup," ...             **\n")
    cat("**                                                   **\n")
    cat("*******************************************************\n\n")
    
    # Create command lines
    CmdLine <- paste("/home/txw/samtools/samtools-1.2/samtools view -h -q1", 
                     BamFileDedup, "-o", BamFileUnique)
    system(CmdLine)
  }
  
  # Run sorting
  if (blnSort){
    cat("\n\n*******************************************************\n")
    cat("**                                                   **\n")
    cat("**    Sorting", BamFileUnique," ...             **\n")
    cat("**                                                   **\n")
    cat("*******************************************************\n\n")
    
    # Command line file for indexing
    CmdLine <- paste("/home/txw/samtools/samtools-1.2/samtools index", 
                     BamFileUnique)
    system(CmdLine)
    
    # Command line file for sorting
    CmdLine <- paste("/home/txw/samtools/samtools-1.2/samtools sort -o", 
                     BamFileSorted2, "-T", paste(BamFile, ".tmp", sep = ""), 
                     BamFileUnique)
    system(CmdLine)
  }
  
  # Create index file for final file
  cat("\n\n*******************************************************\n")
  cat("**                                                   **\n")
  cat("**    Creating index for", BamFileSorted2," ...  **\n")
  cat("**                                                   **\n")
  cat("*******************************************************\n\n")
  
  # Command line file for indexing
  CmdLine <- paste("/home/txw/samtools/samtools-1.2/samtools index", 
                   BamFileSorted2)
  system(CmdLine)
  
}


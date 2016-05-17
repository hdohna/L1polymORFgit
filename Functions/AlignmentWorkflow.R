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

AlignmentWorkflow <- function(InFile, 
                              ReferenceFile, 
                              # SamFile = NULL,
                              # BamFile = NULL,
                              # BamFileSorted = NULL,
                              # BamFileDedup = NULL,
                              # MetricsFileDedup = NULL,
                              # BamFileUnique = NULL,
                              # BamFileSorted2 = NULL,
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
  
  # Function to get and replace file suffix
  GetSuffix <- function(FileName){
    NameSplit <- strsplit(FileName, "\\.")[[1]]
    NameSplit[length(NameSplit)]
  }
  ReplaceSuffix <- function(FileName, Replacement){
    NameSplit <- strsplit(FileName, "\\.")[[1]]
    paste(c(NameSplit[-length(NameSplit)], Replacement), collapse = "")
  }

  # if (is.null(FastqFile) & is.null(SamFile)){
  #   stop("Either FastqFile or SamFile has to be specified!\n")
  # }
  # # Create file names
  # if (is.null(SamFile)) {
  #   FilePrefix <- strsplit(FastqFile, ".fastq")[[1]][1]
  #   SamFile    <- paste(FilePrefix, ".sam", sep = "")
  # }
  # if (is.null(BamFile)) {
  #   FilePrefix <- strsplit(SamFile, ".sam")[[1]][1]
  #   BamFile    <- paste(FilePrefix, ".bam", sep = "")
  # }
  # if (is.null(BamFileSorted)) {
  #   FilePrefix  <- strsplit(BamFile, ".bam")[[1]][1]
  #   BamFileSorted <- paste(FilePrefix, ".sorted.bam", sep = "")
  # }
  # if (is.null(MetricsFileDedup)) {
  #   FilePrefix   <- strsplit(BamFileSorted, ".sorted.bam")[[1]][1]
  #   MetricsFileDedup <- paste(FilePrefix, ".dedup.metrics", sep = "")
  # }
  # if (is.null(BamFileDedup)) {
  #   FilePrefix   <- strsplit(BamFileSorted, ".sorted.bam")[[1]][1]
  #   BamFileDedup <- paste(FilePrefix, ".dedup.bam", sep = "")
  # }
  # if (is.null(BamFileUnique)) {
  #   FilePrefix    <- strsplit(BamFileDedup, ".bam")[[1]][1]
  #   BamFileUnique <- paste(FilePrefix, ".unique.bam", sep = "")
  # }
  # if (is.null(BamFileSorted2)) {
  #   FilePrefix  <- strsplit(BamFileUnique, ".bam")[[1]][1]
  #   BamFileSorted2 <- paste(FilePrefix, ".sorted.bam", sep = "")
  # }
  
  
  # Run alignment
  if (blnAlign){
    cat("\n\n*******************************************************\n")
    cat("**                                                   **\n")
    cat("**    Aligning", InFile," ...             **\n")
    cat("**                                                   **\n")
    cat("*******************************************************\n\n")
    
    if(GetSuffix(InFile) != "fastq"){
      stop("Infile has to be fastq file for alignment!\n")
    }
    OutFile <- ReplaceSuffix(InFile, ".sam")
    CmdLine <- paste(AlignCommand,  ReferenceFile, FastqFile, ">", SamFile)
    system(CmdLine)
    InFile <- OutFile
  }
  
  # Turn sam file into sorted and indexed bam file
  if (blnSam2SortedBam){
    cat("\n\n*******************************************************\n")
    cat("**                                                   **\n")
    cat("**    Converting", InFile," to bam file   **\n")
    cat("**                                                   **\n")
    cat("*******************************************************\n\n")
    
    if(GetSuffix(InFile) != "sam"){
      stop("Infile has to be sam file for conversion to bam file!\n")
    }
    OutFile <- ReplaceSuffix(InFile, ".bam")
    
    # Command line file for conversion
    CmdLine <- paste("/home/txw/samtools/samtools-1.2/samtools view", 
                     InFile, "-b -h -o", OutFile)
    system(CmdLine)
    InFile <- OutFile
    
    cat("\n\n*******************************************************\n")
    cat("**                                                   **\n")
    cat("**    Sorting", InFile,"    **\n")
    cat("**                                                   **\n")
    cat("*******************************************************\n\n")
    
    # Command line file for sorting
    OutFile <- ReplaceSuffix(InFile, ".sorted.bam")
    CmdLine <- paste("/home/txw/samtools/samtools-1.2/samtools sort -o", 
                     OutFile, "-T", paste(InFile, ".tmp", sep = ""), 
                     InFile)
    system(CmdLine)
    InFile <- OutFile
    
    # Command line file for indexing
    CmdLine <- paste("/home/txw/samtools/samtools-1.2/samtools index", 
                     InFile)
    system(CmdLine)
  }

  # Run deduplication
  if (blnDedup){
    cat("\n\n*******************************************************\n")
    cat("**                                                   **\n")
    cat("**    Deduplicating", InFile," ...             **\n")
    cat("**                                                   **\n")
    cat("*******************************************************\n\n")
    
    if(GetSuffix(InFile) != "bam"){
      stop("Infile has to be bam file for deduplication!\n")
    }
    
    # Create command lines
    OutFile          <- ReplaceSuffix(InFile, ".dedup.bam")
    MetricsFile      <- ReplaceSuffix(InFile, ".dedup.metrics")
    DedupInputLine   <- paste('INPUT=', InFile, sep = "")
    DedupOutputLine  <- paste('OUTPUT=', OutFile, sep = "")
    DedupMetricsLine <- paste('METRICS_FILE=', MetricsFile, sep = "")

    CmdLine <- paste(DedupCommand, DedupInputLine, DedupOutputLine, 
                     DedupMetricsLine, DedupOptions)
    system(CmdLine)
    InFile <- OutFile
  }
  
  # Run removal of zero mapping quality
  if (blnRemMapQ0){
    cat("\n\n*******************************************************\n")
    cat("**                                                   **\n")
    cat("**    Remove mapq 0 from", InFile," ...             **\n")
    cat("**                                                   **\n")
    cat("*******************************************************\n\n")
    
    if(GetSuffix(InFile) != "bam"){
      stop("Infile has to be bam file for filtering!\n")
    }
    
    OutFile <- ReplaceSuffix(InFile, ".unique.bam")
    
    # Create command lines
    CmdLine <- paste("/home/txw/samtools/samtools-1.2/samtools view -h -q 1", 
                     InFile, "-o", OutFile)
    system(CmdLine)
    InFile <- OutFile
  }
  
  # Run sorting
  if (blnSort){
    cat("\n\n*******************************************************\n")
    cat("**                                                   **\n")
    cat("**    Sorting", InFile," ...             **\n")
    cat("**                                                   **\n")
    cat("*******************************************************\n\n")
    
    OutFile <- ReplaceSuffix(InFile, ".sorted.bam")

    # Command line file for sorting
    CmdLine <- paste("/home/txw/samtools/samtools-1.2/samtools sort -o", 
                     OutFile, "-T", paste(OutFile, ".tmp", sep = ""), 
                     InFile)
    system(CmdLine)
    InFile <- OutFile
    
  }
  
  
  cat("\n\n*******************************************************\n")
  cat("**                                                   **\n")
  cat("**    Indexing", InFile," ...             **\n")
  cat("**                                                   **\n")
  cat("*******************************************************\n\n")

  # Command line file for indexing
  CmdLine <- paste("/home/txw/samtools/samtools-1.2/samtools index", 
                   InFile)
  system(CmdLine)
  
}


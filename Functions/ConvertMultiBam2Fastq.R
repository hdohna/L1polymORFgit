##############################################
#
# General description:
#
#   The following function aligns all fastq files in a folder
#   to the same reference using bwa

# Input:
#
#     FastQFolder: character string providing path to folder that contains 
#          the fastq files to be mapped
#     Reference: character string providing path to file containing the 
#          common reference
#     IndexCommand: text string providing command to create an index file
#     AlignCommand: text string providing the alignment command (options can be
#          added here)
#     SamSuffix: suffix for sam files created by alignment


# Output:
#   
#     ...

##############################################

ConvertMultiBam2Fastq <- function(InBamFilePaths, NrJobsPerBatch = 100, WaitBetwJobs = 100) {
  
  cat("****************************************************\n")
  cat("**                                                **\n")
  cat("**    Running function ConvertMultiBam2Fastq ...          **\n")
  cat("**                                                **\n")
  cat("****************************************************\n\n")
  
  # Convert into fastq
  OutFastqFilePaths <- gsub(".bam", ".fastq", InBamFilePaths)
  cat("Converting bam to fastq files \n")
  HeaderLines = c('#! /bin/sh','#$ -N TEST', '#$ -cwd',
                  '#$ -j y', '#$ -S /bin/bash', '#', '')
  
  if (length(InBamFilePaths) > NrJobsPerBatch){
    Starts <- seq(1, length(InBamFilePaths), NrJobsPerBatch)
    Ends   <- c(Starts[-1] - 1, length(InBamFilePaths))
  }
  for (j in 1:length(Starts)){
    for (i in Starts[j]:Ends[j]){
      cat("Converting bam file of range", i, "of", length(InBamFilePaths), "\n")
      Bam2FastqCmd <- paste("samtools fastq", InBamFilePaths[i], ">", 
                            OutFastqFilePaths[i]) 
      Commands <- c("module load samtools", Bam2FastqCmd)
      FileName <- paste("/home/hzudohna/tmpBam2Fastq",i, sep = "_")
      ScriptName <- paste("bam2Fastq", i, sep = "_")
      cat("Running commands", paste(Commands, "\n"), "\n")
      CreateAndCallqsubScript(file = FileName, qsubHeaderLines = HeaderLines, 
                              qsubCommandLines = Commands, 
                              scriptName = ScriptName)
    }
    Sys.sleep(WaitBetwJobs)
  }
  
}



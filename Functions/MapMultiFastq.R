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

MapMultiFastq <- function(FastQPaths, 
                          Reference, 
                          IndexCommand,
                          AlignCommand,
                          SamSuffix = "_aln.sam", NrJobsPerBatch = 100, WaitBetwJobs = 100) {
  
  cat("****************************************************\n")
  cat("**                                                **\n")
  cat("**    Running function MapMultiFastq ...          **\n")
  cat("**                                                **\n")
  cat("****************************************************\n\n")
  
  # Get all paths to fastq files in the folder (legacy of previous version)
#  FastQPaths <- list.files(FastQFolder, pattern = ".fastq", full.names = T)
  
  # Create index file
  cat("*******   Creating index for", Reference, "...   *******\n")
  CmdIndex <- c(IndexCommand[1], paste(IndexCommand[2], Reference))
  system(CmdIndex)
  
  # Run BWA for each little fastq file  
  OutFiles <- paste(substr(FastQPaths, 1, nchar(FastQPaths) - 6), SamSuffix, sep = "")
  CmdLines <- paste(AlignCommand[2],  Reference, FastQPaths)
  CmdLines <- paste(CmdLines, OutFiles, sep = " > ")
  HeaderLines = c('#! /bin/sh','#$ -N TEST', '#$ -cwd', '#$ -l h_vmem=3G',
                  '#$ -j y', '#$ -S /bin/bash', '#', '')

  if (length(CmdLines) > NrJobsPerBatch){
    Starts <- seq(1, length(CmdLines), NrJobsPerBatch)
    Ends   <- c(Starts[-1] - 1, length(CmdLines))
  } else {
    Starts <- 1
    Ends <- length(CmdLines)
  }
  for (j in 1:length(Starts)){
    for (i in Starts[j]:Ends[j]){
      CmdLocal <- c(AlignCommand[1], CmdLines[i])
    FileName <- paste("/home/hzudohna/tmpBWA",i, sep = "_")
    ScriptName <- paste("bwa", i, sep = "_")
    cat("Running commands\n", paste(CmdLocal, "\n"), "\n")
    # system(CmdLocal)
    CreateAndCallqsubScript(file = FileName, qsubHeaderLines = HeaderLines,
                            qsubCommandLines = CmdLocal,
                            scriptName = ScriptName)
    
    }

    Sys.sleep(WaitBetwJobs)
  }
  # Return paths to fastq files and sam files
  cbind.data.frame(FastQPath = FastQPaths, SamPath = OutFiles)
}



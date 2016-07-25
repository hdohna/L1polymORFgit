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

MapMultiFastq <- function(FastQFolder, Reference, 
                          IndexCommand,
                          AlignCommand,
                          SamSuffix = "_aln.sam") {
  
  cat("****************************************************\n")
  cat("**                                                **\n")
  cat("**    Running function MapMultiFastq ...          **\n")
  cat("**                                                **\n")
  cat("****************************************************\n\n")
  
  # Get all paths to fastq files in the folder
  FastQPaths <- list.files(FastQFolder, pattern = ".fastq", full.names = T)
  
  # Create index file
  cat("*******   Creating index for", Reference, "...   *******\n")
  CmdIndex <- c(IndexCommand[1], paste(IndexCommand[2], Reference))
  browser()
  system(CmdIndex)
  
  # Run BWA for each little fastq file  
  OutFiles <- paste(substr(FastQPaths, 1, nchar(FastQPaths) - 6), SamSuffix, sep = "")
  CmdLines <- paste(AlignCommand[2],  Reference, FastQPaths)
  CmdLines <- c(AlignCommand[1], paste(CmdLines, OutFiles, sep = " > "))
  for (CmdL in CmdLines) {
    cat("Running command", CmdL, "\n")
    system(c(AlignCommand[1], CmdL))
  }
  
  # Return paths to fastq files and sam files
  cbind.data.frame(FastQPath = FastQPaths, SamPath = OutFiles)
}



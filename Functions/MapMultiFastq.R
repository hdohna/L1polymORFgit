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
                          IndexCommand = '/home/txw/bwa/bwa-0.7.12/bwa index',
                          AlignCommand = '/home/txw/bwa/bwa-0.7.12/bwa mem',
                              SamSuffix = "_aln.sam") {
  
  # Bam suffices for alignment files created 
  BamSuffix <- paste(substr(SamSuffix, 1, nchar(SamSuffix) - 4), ".bam", sep = "")
  
  cat("*******   Mapping little fastq files in", FastQFolder, "...   *******\n")
  
  # Get all paths to fastq files in the folder
  FastQPaths <- list.files(FastQFolder, pattern = ".fastq", full.names = T)
  
  # Create index file
  CmdIndex <- paste(IndexCommand, Reference)
  system(CmdIndex)
  
  # Run BWA for each little fastq file  
  OutFiles <- paste(substr(FastQPaths, 1, nchar(FastQPaths) - 6), SamSuffix, sep = "")
  CmdLines <- paste(AlignCommand,  Reference, FastQPaths)
  CmdLines <- paste(CmdLines, OutFiles, sep = " > ")
  for (CmdL in CmdLines) system(CmdL)
  
  # Return paths to fastq files and sam files
  cbind.data.frame(FastQPath = FastQPaths, SamPath = OutFiles)
}



##############################################
#
# General description:
#
#   The following function adds read groups to all bam files in a folder
#   

# Input:
#
#     BamFolder: character string providing path to folder that contains 
#          the fastq files to be mapped
#     Reference: character string providing path to file containing the 
#          common reference
#     AddGroupCmd: text string providing command to add read groups
#     AlignCommand: text string providing the alignment command (options can be
#          added here)
#     ReadGroupSuffix: suffix for bam files created after adding read groups


# Output:
#   
#     ...

##############################################

CreateMultiBamIndex <- function(BamFolder,  
   CreateBamIndexCmd   = "/home/txw/samtools/samtools-1.2/samtools index",
   BamSuffix = ".bam"
) {
  
  cat("*******************************************************\n")
  cat("**                                                   **\n")
  cat("**    Running function CreateMultiBamIndex ...       **\n")
  cat("**                                                   **\n")
  cat("*******************************************************\n\n")
  
  # Get all paths to sam files in the folder
  FileNames <- list.files(BamFolder, pattern = BamSuffix, full.names = T)
  Bam.Pattern <- grep(".bam.", FileNames)
  if (length(Bam.Pattern) > 0){
    FileNames <- FileNames[-grep(".bam.", FileNames)]
  }
  
  # Create a command per file 
  CmdLines <- paste(CreateBamIndexCmd,  FileNames)
  
  # Loop over command line and run them 
  for (CmdL in CmdLines) {
    cat("Calling command:\n", CmdL, "\n")
    system(CmdL)
  }
  
}



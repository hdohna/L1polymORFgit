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

AddMultiReadGroups <- function(BamFolder,  
   AddGroupCmd   = "java -jar /home/txw/picard/picard-tools-1.131/picard.jar AddOrReplaceReadGroups",
   AddGroupOptions = c("RGLB=lib1", "RGPL=illumina", "RGPU=unit1", "RGSM=20",
                    "SORT_ORDER=null", "CREATE_INDEX=TRUE", "VALIDATION_STRINGENCY=LENIENT"),
   ReadGroupSuffix = "withRG.bam") {
  
  cat("*******************************************************\n")
  cat("**                                                   **\n")
  cat("**    Running function AddMultiReadGroups ...        **\n")
  cat("**                                                   **\n")
  cat("*******************************************************\n\n")
  
  # Get all paths to fastq files in the folder
  FileNames <- list.files(BamFolder, pattern = ".bam", full.names = T)
  FileNames <- FileNames[-grep(".bam.", FileNames)]
  FileNames <- FileNames[-grep(ReadGroupSuffix, FileNames)]
  
  # Create output file names
  OutFileNames <- substr(FileNames, 1, nchar(FileNames) - 4)
  OutFileNames <- paste(OutFileNames, ReadGroupSuffix, sep = "_")
  
  # Create a command per file 
  OptionLines <- paste(AddGroupOptions, collapse = " ")
  InFiles  <- paste("I=", FileNames, sep = "")
  OutFiles <- paste("O=", OutFileNames, sep = "")
  CmdLines <- paste(AddGroupCmd,  InFiles, OutFiles, OptionLines)
  
  # Loop over command line and run them 
  for (CmdL in CmdLines) system(CmdL)
  
  # Return paths to fastq files and sam files
  cbind.data.frame(InFileNames = FileNames, OutFileNames = OutFileNames)
}



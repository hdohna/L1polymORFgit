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
#     IndexCommand: text string providing command to create an index file
#     AlignCommand: text string providing the alignment command (options can be
#          added here)
#     SamSuffix: suffix for sam files created by alignment


# Output:
#   
#     ...

##############################################

CallMultiVariants <- function(BamFolder,  
 HapTypeCallCmd = "java -jar /home/txw/GATK/GenomeAnalysisTK-2.1-11-g13c0244/GenomeAnalysisTK.jar -T HaplotypeCaller",
    RefSeqPath,
    OptionLines = "--emitRefConfidence GVCF",
    BamSuffix = "withRG.bam",
    VCFSuffix = ".vcf") {
    
  cat("*******   Calling haplotypes for files in", BamFolder, "...   *******\n")
  
  # Get all paths to fastq files in the folder
  FileNames <- list.files(BamFolder, pattern = BamSuffix, full.names = T)
  Bam.Pattern <- grep(".bam.", FileNames)
  if (length(Bam.Pattern) > 0){
    FileNames <- FileNames[-grep(".bam.", FileNames)]
  }

  # Create output file names
  OutFileNames <- substr(FileNames, 1, nchar(FileNames) - 4)
  OutFileNames <- paste(OutFileNames, VCFSuffix, sep = "")
  
  # Create a command per file 
  OptionLines <- paste(OptionLines, collapse = " ")
  RefFile  <- paste("-R", RefSeqPath)
  InFiles  <- paste("-I", FileNames)
  OutFiles <- paste("-o", OutFileNames)
  CmdLines <- paste(HapTypeCallCmd,  RefFile, InFiles, OutFiles, OptionLines)
  
  # Loop over command line and run them 
  for (CmdL in CmdLines) {
    cat("Calling command:\n", CmdL, "\n")
    system(CmdL)
    }
  
  # Return paths to fastq files and sam files
  cbind.data.frame(InFileNames = FileNames, OutFileNames = OutFileNames)
}



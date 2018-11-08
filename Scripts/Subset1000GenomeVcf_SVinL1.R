# The script below identifies deletions that are likely L1 polymorphis,s

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specify parameters
DataFolder     <- "/labs/dflev/hzudohna/1000Genomes/"
LineIdentifier <- "SVTYPE=DEL"
FileSuffix     <- "ME_Deletions.vcf"

# Get file names, loop over files and do the filtering
# Example file name: ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
AllFiles <- list.files(DataFolder, pattern = "phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf", 
                       full.names = T)
AllFiles <- AllFiles[-grep("vcf.", AllFiles)]
InFile <- AllFiles[1]
for (InFile in AllFiles){
  InFileSplit <- strsplit(InFile, "\\.")[[1]]
  OutFile     <- paste(paste(InFileSplit[1:2], collapse = "_"), FileSuffix, sep = "_")
  BedLoadCmd  <- "module load bedtools" 
  SubsetCmd   <- paste("bedtools intersect -a", InFile,
                     "-b /labs/dflev/hzudohna/1000Genomes/L1HSRefRanges_Plus200_hg19.bed -wa | ",
                     "grep", LineIdentifier,  ">", OutFile)
  cat("Filtering", InFile, "\n\n")
  ScriptName <- paste("ME_DEL_Script", InFileSplit[2], sep = "_")
  CreateAndCallSlurmScript(file = ScriptName, SlurmCommandLines = c(BedLoadCmd, SubsetCmd))
}

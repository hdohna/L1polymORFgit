# The script below identifies deletions that are likely L1 polymorphis,s

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specify parameters
DataFolder     <- "/labs/dflev/hzudohna/1000Genomes/"
LineIdentifier <- "SVTYPE=DEL"
OutPrefix     <- "L1NoPoly_Indels"

# Get file names, loop over files and do the filtering
# Example file name: ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
AllFiles <- list.files(DataFolder, pattern = "phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf", 
                       full.names = T)
AllFiles <- AllFiles[-grep("vcf.", AllFiles)]
InFile <- AllFiles[1]
for (InFile in AllFiles){
  InFileSplit <- strsplit(InFile, "\\.")[[1]]
  OutPath     <- paste(c(InFileSplit[1:2], OutPrefix), collapse = "_")
  OutPath2    <- paste(c(InFileSplit[2], OutPrefix), collapse = "_")
  VcfLoadCmd  <- "module load vcftools" 
  SubsetCmd   <- paste("vcftools --vcf", InFile,
                     "--bed /labs/dflev/hzudohna/RefSeqData/L1GRangesNoPoly.bed --keep-only-indels",
                     "--recode --recode-INFO-all",
                     "--stdout >", OutPath)
  cat("Filtering", InFileSplit[2], "\n\n")
  ScriptName <- paste("L1indelL_Script", InFileSplit[2], sep = "_")
  CreateAndCallSlurmScript(file = ScriptName, SlurmCommandLines = c(VcfLoadCmd, SubsetCmd))
}

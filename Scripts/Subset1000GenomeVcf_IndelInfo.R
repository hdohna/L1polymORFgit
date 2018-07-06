# The script below subsets vcf files from the 1000 genome project to retain
# only rows of indels and the info columns

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specify parameters
DataFolder     <- "/labs/dflev/hzudohna/1000Genomes/"
LineIdentifier <- "VT=INDEL"
FilePrefix     <- "IndelInfo"
LineIdentifier <- "ME:LINE1"
FilePrefix     <- "LINE1_Info"

# Get file names, loop over files and do the filtering
# Example file name: ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
AllFiles <- list.files(DataFolder, pattern = "genotypes.vcf", full.names = T)
AllFiles <- AllFiles[-grep("vcf.", AllFiles)]
InFile <- AllFiles[1]
for (InFile in AllFiles){
  InFileSplit <- strsplit(InFile, "\\.")[[1]]
  OutFile     <- paste(c(InFileSplit[1], FilePrefix, "_", InFileSplit[2], ".vcf"), collapse = "")
  OutFile     <- gsub("ALL", "", OutFile)
  OutFile
  SubsetCmd <- paste("cut -s -f 1-9 ", InFile, 
                     paste("| grep", LineIdentifier,  ">", OutFile))
  cat("Filtering", InFile, "\n\n")
  ScriptName <- paste("IndelScript", InFileSplit[2], sep = "")
  CreateAndCallSlurmScript(file = ScriptName, SlurmCommandLines = SubsetCmd)
}

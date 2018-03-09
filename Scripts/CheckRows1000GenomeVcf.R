# The script below checks rows in vcf files from the 1000 genome project to 
# identify rows with unusual numbers of characters
# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specify parameters
DataFolder     <- "/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/"
FilePrefix     <- "weirdLines"

# Get file names, loop over files and do the filtering
# Example file name: ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
AllFiles <- list.files(DataFolder, pattern = "ALL.chr", full.names = T)
AllFiles <- AllFiles[grep(".vcf", AllFiles)]
AllFiles <- AllFiles[-grep("vcf.", AllFiles)]
InFile <- AllFiles[1]
for (InFile in AllFiles){
  InFileSplit <- strsplit(InFile, "\\.")[[1]]
  OutFile     <- paste(InFileSplit[1:2], collapse = "_")
  OutFile     <- gsub("ALL", FilePrefix, OutFile)
  Cmd <- paste("qsub pbs_checkLLVcf1000G", InFile,  OutFile)
  cat("Checking lines in", InFile, "\n\n")
  system(Cmd)
}

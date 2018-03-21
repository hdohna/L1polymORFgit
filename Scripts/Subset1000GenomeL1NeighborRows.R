# The script below subsets rows neighboring L1 entries in
# 1000 genome vcf files

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specify parameters
DataFolder     <- "/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/"
NrRows         <- 100
FilePrefix     <- "NeighborLinesL1"
FilePattern    <- "ALL"
ScriptName     <- "pbs_AdjLines" 

# Get file names, loop over files and do the filtering
# Example file name: ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
AllFiles <- list.files(DataFolder, pattern = FilePattern, full.names = T)
AllFiles <- AllFiles[-grep("vcf.", AllFiles)]
AllFiles <- AllFiles[-grep(".ped", AllFiles)]
AllFiles <- AllFiles[-grep("panel", AllFiles)]
AllFiles 
InFile <- AllFiles[1]
InFile
for (i in 2:length(AllFiles)){
  InFile <- AllFiles[i]
  InFileSplit <- strsplit(InFile, "\\.")[[1]]
  OutFile     <- paste(c(InFileSplit[1:2], ".vcf"), collapse = "")
  OutFile     <- gsub("ALL", FilePrefix, OutFile)
  Cmd <- paste("qsub", ScriptName,NrRows, InFile, OutFile)
  cat("Filtering", InFile, "\n\n")
  system(Cmd, wait = F)
}

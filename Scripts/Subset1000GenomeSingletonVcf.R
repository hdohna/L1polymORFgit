# The script below subsets vcf files from the 1000 genome project to retain
# only files with singletons (variant that occurred only once)

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specify parameters
DataFolder     <- "/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/"
FilePrefix     <- "Singleton"

# Get file names, loop over files and do the filtering
AllFiles <- list.files(DataFolder, pattern = ".vcf", full.names = T)
AllFiles <- AllFiles[-grep("vcf.", AllFiles)]
AllFiles <- AllFiles[grep("ALL.", AllFiles)]
InFile <- AllFiles[11]
for (InFile in AllFiles){
  InFileSplit <- strsplit(InFile, "\\.")[[1]]
  OutFile     <- paste(InFileSplit[1:2], collapse = "_")
  OutFile     <- gsub("ALL", FilePrefix, OutFile)
  Cmd <- paste("qsub pbs_Singleton", InFile, OutFile)
  cat("Filtering", InFile,"\n")
  cat("Results will be saved to", OutFile,"\n\n")
  system(Cmd, wait = F)
}

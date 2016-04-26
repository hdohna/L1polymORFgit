# The script below subsets vcf files from the 1000 genome project to retain
# only files with Line1 insertions

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Specify parameters
DataFolder     <- "/share/diskarray3/hzudohna/1000Genomes/"
LineIdentifier <- "INS:ME:LINE1"
FilePrefix     <- "LINE1"

# Get file names, loop over files and do the filtering
AllFiles <- list.files(DataFolder, pattern = ".vcf", full.names = T)
InFile <- AllFiles[1]
for (InFile in AllFiles){
  InFileSplit <- strsplit(InFile, "\\.")[[1]]
  OutFile     <- paste(c(InFileSplit[1:2], ".vcf"), collapse = "")
  OutFile     <- gsub("ALL", FilePrefix, OutFile)
  Cmd <- paste("grep", LineIdentifier,  InFile, ">", OutFile)
  cat("Filtering", InFile, "\n\n")
  system(Cmd)
}

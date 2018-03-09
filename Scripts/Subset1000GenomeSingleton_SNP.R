# The script below subsets singleton files from the 1000 genome project to retain
# only files with SNPs

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specify parameters
DataFolder     <- "/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/"
LineIdentifier <- "INS:ME:LINE1"
FilePrefix     <- "L1"
LineIdentifier <- "VT=SNP"
FilePrefix     <- "SNP"

# Get file names, loop over files and do the filtering
# Example file name: ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
AllFiles <- list.files(DataFolder, pattern = "Singleton_chr", full.names = T)
InFile <- AllFiles[1]
InFile
for (i in 1:length(AllFiles)){
  InFile <- AllFiles[i]
  InFileSplit <- strsplit(InFile, "\\_")[[1]]
  OutFile     <- paste(c(InFileSplit[1], FilePrefix, InFileSplit[2]), collapse = "_")
  OutFile     <- gsub("ALL", FilePrefix, OutFile)
  Cmd <- paste("grep", LineIdentifier,  InFile, ">", OutFile)
  cat("Filtering", InFile, "\n\n")
  CreateAndCallqsubScript_scg4(file = paste("filterSNP", i, sep = ""),
     qsubCommandLines = Cmd, scriptName = 'filterSNP')
}

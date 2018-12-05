# The script below vcf files from the 1000 Genome data and identifies columns containing
# unphased genotypes
# Data were downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specify parameters
NrInfoCols <- 9
DataFolder <- "/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/"
OutFile    <- paste(DataFolder, "ColnamesUnphased", sep = "")

# Get column names from vcf files
cat("Adding column names to concatenated L1 table\n")
VcfFiles <- list.files(DataFolder, pattern = "genotypes.vcf", full.names = T)
VcfFiles <- VcfFiles[-grep("vcf.", VcfFiles)]
MEI1000GLines <- readLines(VcfFiles[1], n = 300)
StartChar <- substr(MEI1000GLines, 1, 1)
StartLines <- MEI1000GLines[StartChar == "#"]
FileLines  <- apply(L1Table, 1, function(x) paste(x, collapse = "\t"))

# Add column names
ColNames  <- MEI1000GLines[max(which(StartChar == "#"))]
ColNames  <- gsub("#", "", ColNames)
ColNames  <- strsplit(ColNames, "\t")[[1]]

# Loop over file names, read file and append to existing
cat("Reading vcf files by chromosome\n")
VcfFile <- VcfFiles[1]
Cols2Remove <- c()
for (VcfFile in VcfFiles){
  cat("Processing", VcfFile, "\n")
  NewTable <- read.delim(VcfFile, skip = length(StartLines),
                         sep = "\t",header = F, as.is = T)
  NewTable[49, 1:12]
  blnCol2 <- apply(NewTable[-1,], 2, FUN = function(x) 
    length(grep("2", as.character(x))) > 0)
  blnRow2 <- apply(NewTable[,10:ncol(NewTable)], 1, FUN = function(x) 
    length(grep("2", as.character(x))) > 0)
  idxCol2 <- which(blnCol2)
  idxCol2 <- idxCol2[idxCol2 > 9]
  idxCol2
  idxRow2 <- which(blnRow2)
  idxRow2
  NewTable[49, 9 + grep("2", NewTable[49, 10:ncol(NewTable)])]
  nrow(NewTable)
  length(idxRow2)
}


# Write out concatenated L1 file
writeLines(c(ColNames, FileLines), OutFileAllL1)

# Add column names
ColNames  <- MEI1000GLines[max(which(StartChar == "#"))]
ColNames  <- strsplit(ColNames, "\t")[[1]]
ColNames  <- gsub("#", "", ColNames)
colnames(L1Table) <- ColNames

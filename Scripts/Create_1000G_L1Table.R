# The script below subsets creates a table with L1s from the 1000 Genome data
# Data were downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specify parameters
NrInfoCols <- 9
DataFolder <- "/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/"
FilePrefix <- "LINE1"
OutFileAllL1  <- paste(DataFolder, "L1all.vcf", sep = "")
OutFile       <- paste(DataFolder, "L1_1000G_withGenoNum", sep = "")

# function to get the genotype
GetGenotype <- function(SampleColumn){
  sapply(as.character(SampleColumn), function(x){
    Split <- strsplit(x, "\\|")[[1]]
    sum(as.numeric(Split))
  })
}

# Get file names, loop over files and do the filtering
AllFiles <- list.files(DataFolder, pattern = FilePrefix, full.names = T)
AllFiles <- AllFiles[-grep("chL", AllFiles)]
AllFiles <- AllFiles[-grep("all", AllFiles)]
AllFiles

# Loop over file names, read file and append to existing
cat("Reading vcf files by chromosome\n")
L1Table <- read.delim(AllFiles[1], header = F, skip = 1)
for (L1File in AllFiles[-1]){
  cat("Processing", L1File, "\n")
  NewTable <- read.delim(L1File, header = F, skip = 1)
  L1Table  <- rbind(L1Table, NewTable)
}

# Get column names from vcf files
cat("Adding column names to concatenated L1 table\n")
VcfFiles <- list.files(DataFolder, pattern = "genotypes.vcf", full.names = T)
VcfFiles <- VcfFiles[-grep("vcf.", VcfFiles)]
MEI1000GLines <- readLines(VcfFiles[1], n = 300)
StartChar <- substr(MEI1000GLines, 1, 1)
StartLines <- MEI1000GLines[StartChar == "#"]
FileLines  <- apply(L1Table, 1, function(x) paste(x, collapse = "\t"))

# Write out concatenated L1 file
writeLines(c(StartLines, FileLines), OutFileAllL1)

# Add column names
ColNames  <- MEI1000GLines[max(which(StartChar == "#"))]
ColNames  <- strsplit(ColNames, "\t")[[1]]
ColNames  <- gsub("#", "", ColNames)
colnames(L1Table) <- ColNames

# Add column names and turn genotypes into numeric values
# function to get the genotype
GetGenotype <- function(SampleColumn){
  sapply(as.character(SampleColumn), function(x){
    Split <- strsplit(x, "\\|")[[1]]
    sum(as.numeric(Split))
  })
}

# Turn genotype columns into numeric values
cat("Converting genotypes into numeric values\n")
SampleColumns <- colnames(L1Table)[(NrInfoCols + 1):ncol(L1Table)]
GenoDF <- sapply(SampleColumns, function(x) GetGenotype(L1Table[,x]))
colnames(GenoDF) <- SampleColumns
L1Table <- cbind(L1Table[,1:NrInfoCols], GenoDF)

# Extract insertion length
cat("Extracting insertion length\n")
L1Table$InsLength <- sapply(as.character(L1Table$INFO), function(x){
  Xsplit <- strsplit(x, ";")[[1]]
  SVLENchar <- grep("SVLEN", Xsplit, value = T)
  if (length(SVLENchar) > 0){
    as.numeric(substr(SVLENchar, 7, nchar(SVLENchar)))
  } else {
    NA
  }
})

# Write out L1Table
cat("Writing out table\n")
write.table(L1Table, OutFile)




# The script below subsets creates a table with L1s from the 1000 Genome data

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specify parameters
NrInfoCols <- 9
DataFolder <- "/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/"
FilePrefix <- "LINE1"
OutFile    <- paste(DataFolder, "L1_1000G_withGeno")

# function to get the genotype
GetGenotype <- function(SampleColumn){
  sapply(as.character(SampleColumn), function(x){
    Split <- strsplit(x, "\\|")[[1]]
    sum(as.numeric(Split))
  })
}

# Get file names, loop over files and do the filtering
AllFiles <- list.files(DataFolder, pattern = FilePrefix, full.names = T)

# Loop over file names, read file and append to existing
L1Table <- read.delim(AllFiles[1], header = F)
for (L1File in AllFiles[-1]){
  NewTable <- read.delim(L1File, header = F)
  L1Table  <- rbind(L1Table, NewTable)
}

# Get column names from vcf files
VcfFiles <- list.files(DataFolder, pattern = "genotypes.vcf", full.names = T)
VcfFiles <- VcfFiles[-grep("vcf.", VcfFiles)]
MEI1000GLines <- readLines(VcfFiles[1], n = 300)
StartChar <- substr(MEI1000GLines, 1, 2)
ColNames  <- MEI1000GLines[max(which(StartChar == "##")) + 1]
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
SampleColumns <- colnames(L1Table)[(NrInfoCols + 1):ncol(L1Table)]
GenoDF <- sapply(SampleColumns, function(x) GetGenotype(L1Table[,x]))
colnames(GenoDF) <- SampleColumns
L1Table <- cbind(L1Table[,1:NrInfoCols], GenoDF)

# Extract insertion length
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
write.table(L1Table, OutFile)




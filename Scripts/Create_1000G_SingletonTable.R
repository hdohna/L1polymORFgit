# The script below subsets creates a table with L1s from the 1000 Genome data
# Data were downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')
library(data.table)

# Specify parameters
NrInfoCols  <- 9
NrSamples   <- 2504
TotCols     <- NrInfoCols + NrSamples
NLines2Read <- 10^3
DataFolder <- "/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/"
FilePrefix <- "Singleton"

# Get column names from vcf files
cat("Adding column names to concatenated L1 table\n")
VcfFiles <- list.files(DataFolder, pattern = "genotypes.vcf", full.names = T)
VcfFiles <- VcfFiles[-grep("vcf.", VcfFiles)]

File       <- VcfFiles[1]
StartLines <- readLines(File, n = 300)
StartChar  <- substr(StartLines, 1, 2)
nSkip      <- sum(StartChar == "##") + 1
ColNames   <- StartLines[nSkip]
ColNames   <- strsplit(ColNames, "\t")[[1]]
ColNames   <- gsub("#", "", ColNames)
SampleNames <- ColNames[NrInfoCols + 1:NrSamples]

Nread      <- NLines2Read
FileSplit <- strsplit(File, "\\.")[[1]]
OutFile <- paste(FilePrefix, FileSplit[2], FileSplit[length(FileSplit)],
                 sep = ".")
OutPath    <- paste(DataFolder, OutFile, sep = "")
NewCon     <- file(OutPath)
TotRead    <- 0
TotWritten <- 0

while(Nread == NLines2Read){
  CurrentLines <- scan(File, skip = nSkip, nlines = NLines2Read,
                       sep = "\t", what = character())
  idxNewLines <- seq(1, length(CurrentLines), TotCols)
  Nread       <- length(idxNewLines)
  TotRead <- TotRead + Nread
  cat("Total lines read:", TotRead, "\n")
  idxRemove   <- rep(idxNewLines, each = NrInfoCols) + 
    rep(0:(NrInfoCols - 1), length(idxNewLines))
  CurrentGeno <- CurrentLines[-idxRemove]
  BarSplit    <- strsplit(CurrentGeno, split = "\\|")
  BarSplitN   <- as.numeric(unlist(BarSplit))
  NGenoMat    <- matrix(BarSplitN, nrow = 2*NrSamples)
  blnSingle   <- colSums(NGenoMat) == 1
  NewLines    <- sapply(which(blnSingle), function(x){
    NL <- c(CurrentLines[idxNewLines[x] + 0:(NrInfoCols - 1)],
            SampleNames[NGenoMat[,x] == 1])
    paste(paste(NL, collapse = "\t"), "\n", sep = "") 
  })
  open(NewCon, open = "w")
  writeLines(NewLines, con = NewCon)
  TotWritten <- TotWritten + length(NewLines)
  cat("Total lines written:", TotWritten, "\n")
  
}
close(NewCon)
cat("All singleton written out\n")

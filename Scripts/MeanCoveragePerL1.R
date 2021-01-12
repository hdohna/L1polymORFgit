# The following script calculates the average 1000 genome coverage for each bp 
# within reference L1

# Source start script
source('/home/hb54/L1polymORFgit/Scripts/_Start_L1polymORF_AUB.R')

# Path to reference data and for 1000 genomes
#DataFolder         <- "/labs/dflev/hzudohna/1000Genomes/"
DataFolder         <- '/home/hb54/1000GenomeData/'
#load('/labs/dflev/hzudohna/1000Genomes/GRanges_L1_1000Genomes.RData')

# Get names of potential coverage files
#CoverFiles <- paste(DataFolder, "L1Coverage_", SampleColumns, sep = "")
CoverFiles <- list.files(DataFolder, pattern = "L1Coverage_", full.names = T)

# Loop over files that contain L1 coverage of individual genomes and sum
# coverage values and squared coverage values
NrFiles2Analyze <- 100
CoverSum   <- 0
CoverSumSq <- 0
SuccessiveCors <- NULL
NLines <- NULL
#for (CoverFile in CoverFiles[file.exists(CoverFiles)][1:NrFiles2Analyze]){
for (CoverFile in CoverFiles){
  cat("*******   Getting coverage for", CoverFile, "    *******\n")
  CoverLines <- readLines(CoverFile)
  NLines <- c(length(CoverLines), NLines)
  cat("Number of lines in file:", length(CoverLines), "\n")
  if (length(CoverLines) == 284260){
    NewCover   <- read.table(text = CoverLines)[,3]
    if (CoverFile != CoverFiles[1]){
      SuccessiveCors <- c(SuccessiveCors, cor(OldCover, NewCover))
    }
    OldCover   <- NewCover
    CoverSum   <- CoverSum   + NewCover
    CoverSumSq <- CoverSumSq + NewCover^2
  }
}

# Create a file with mean coverage and mean squared coverage
L1CoverTable <- read.table(CoverFile)[,1:2]
colnames(L1CoverTable) <- c("Chromosome", "Pos")
L1CoverTable$CoverMean <- CoverSum / sum(NLines == 284260)
L1CoverTable$CoverVar <- (CoverSumSq / sum(NLines == 284260)) - L1CoverTable$CoverMean^2

# Write results 
OutPath <- paste(DataFolder, "L1CoverageResults_HG38.RData", sep = "")
cat("Saving results to", OutPath, "   ....")
save.image(file = OutPath)
cat("done!\n")

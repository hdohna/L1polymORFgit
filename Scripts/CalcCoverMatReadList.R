# The script below reads in reads aligned to L1, creates a list of reads, 
# calculates a coverage and quantile matrix

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(ShortRead)
library(Rsamtools)
library(rtracklayer)
library(csaw)
library(GenomicRanges)

# Set file paths
CoveragePlotPath     <- 'D:/L1polymORF/Figures/L1InsertionCoverage_NA12878_PacBio_New.pdf'
OutputFilePath       <- 'D:/L1polymORF/Data/L1_NA12878_PacBio_Coverage_New.RData'
OutFolderName_NonRef <- "D:/L1polymORF/Data/BZ_NonRef"

# Load ranges
load("D:/L1polymORF/Data/BZ_L1Ranges.RData")

# get names of newly created bam files
FileNames <- list.files(OutFolderName_NonRef, pattern = ".bam",
                        full.names = T)
FileNames <- FileNames[-grep(".bam.", FileNames)]

# Loop through file names and read in bam files of reads mapped to L1
cat("*******  Scanning bam files per L1   *************\n")
ScannedL1Ranges <- lapply(1:length(FileNames), function(i) {
  cat("Scanning bam file", i, "of", length(FileNames), "\n")
  scanBam(FileNames[i])
})

# Count the number of reads mapped
NrMapped2L1 <- sapply(ScannedL1Ranges, function(x){
  sum(!is.na(x[[1]]$pos))
})
idxFilesWithReads <- which(NrMapped2L1 > 0)
FilesWithReads    <- FileNames[idxFilesWithReads]

# Get read list per peak
ReadListPerPeak <- lapply(idxFilesWithReads, function(x) {
  ScannedL1Ranges[[x]][[1]]
})
cat("*******  Calculating coverage matrix   *************\n")
CoverMat <- t(sapply(1:length(ReadListPerPeak), function(i) {
  cat("Analyzing row", i, "of", length(ReadListPerPeak), "\n")
  x <- ReadListPerPeak[[i]]
  primMap <- x$flag <= 2047
  RL <- lapply(x, function(y) y[primMap])
  CoverageFromReadList(RL, End = 6064)
}))
dim(CoverMat)
rownames(CoverMat) <- FilesWithReads

# Plot mean coverage
plot(colMeans(CoverMat), ylab = "Mean coverage", xlab = "Position on L1", type = "s")
CreateDisplayPdf('D:/L1polymORF/Figures/L1InsertionCoverage_NA12878_PacBio.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

# Get means and 95% quantiles 
cat("*******  Calculating quantile matrix   *************\n")
QuantileMat <- apply(CoverMat, 2, FUN = function(x) quantile(x, c(0.05, 0.5, 0.95)))
idxFw <- 1:ncol(CoverMat)
idxRv <- ncol(CoverMat):1
plot(QuantileMat[2,], type = "n", ylim = c(0, max(QuantileMat)), 
     ylab = 'Coverage', xlab = "Genomic position")
polygon(c(idxFw, idxRv), c(QuantileMat[1, idxFw], QuantileMat[3, idxRv]),
        col = "grey", border = NA)
lines(QuantileMat[2,], lwd = 1.2)

# Save objects to file
cat("*****  Saving data to", OutputFilePath,  "*********\n")
save(list = c("ScannedL1Ranges", "idxFilesWithReads", "ReadListPerPeak", 
              "CoverMat", "QuantileMat", "FilesWithReads"), file = OutputFilePath)
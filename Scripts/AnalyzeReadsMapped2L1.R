# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(ShortRead)
library(Rsamtools)
library(rtracklayer)
library(csaw)
library(GenomicRanges)

# Load ranges
load("D:/L1polymORF/Data/BZ_L1Ranges.RData")

# Read in table with known L1 
L1Catalogue <- read.csv("D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv", as.is = T)

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
L1CatalogL1Mapped <- L1Catalogue[which(blnL1Mapped & blnAllele1),]

# Lift coordinates and get genomic ranges for catalog L1 on hg19
LiftOverList <- LiftoverL1Catalog(L1CatalogL1Mapped,
                                  ChainFilePath = "D:/L1polymORF/Data/hg38ToHg19.over.chain")
L1CatalogGR <- LiftOverList$GRCatalogue_hg19# Specify folder 
OutFolderName_NonRef <- "D:/L1polymORF/Data/BZ_NonRef"

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
FilesWithReads <- FileNames[idxFilesWithReads]

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


# # Determine range index from file name 
plot(CoverMat[1,], type = "s", xlab = "Position on L1",
     ylab = "Coverage", ylim = c(0, 100))
Cols <- rainbow(nrow(CoverMat))
for (i in 1:nrow(CoverMat)){
  lines(CoverMat[i,], type = "s", col = Cols[i])
}

# Collect information on insertion that fullfill a certain minimum criterion
idx5P <- which(sapply(1:nrow(CoverMat), function(x) all(CoverMat[x, 50:2000] > 0)))
idxFull <- which(sapply(1:nrow(CoverMat), function(x) all(CoverMat[x, 50:6040] > 0)))
FullL1Info <- t(sapply(FilesWithReads[idx5P], function(x){
  FPathSplit <- strsplit(x, "/")[[1]]
  FName      <- FPathSplit[length(FPathSplit)]
  FName      <- substr(FName, 1, nchar(FName) - 4)
  strsplit(FName, "_")[[1]]
}))
FullL1Info <- base::as.data.frame(FullL1Info, stringsAsFactors = F)
colnames(FullL1Info) <- c("chromosome", "idx")
FullL1Info$Max3P <- sapply(idx5P, function(x) max(which(CoverMat[x, ] > 0)))
FullL1Info$idx   <- as.numeric(FullL1Info$idx)
FullL1Info$start <- start(IslGRanges_reduced[FullL1Info$idx])
FullL1Info$end   <- end(IslGRanges_reduced[FullL1Info$idx])
FullL1Info$cover   <- CoverMat[idx5P, 100]
FilesWithReads[idxFull]

# Create genomic ranges for L1 insertions
GRL1Capture <- makeGRangesFromDataFrame(FullL1Info)

load("D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData")

sum(overlapsAny(GRL1Capture, L1CatalogGR))
sum(overlapsAny(GRL1Capture, GRL1Ins1000G))
sum(overlapsAny(IslGRanges_reduced, L1CatalogGR))
sum(overlapsAny(IslGRanges_reduced, GRL1Ins1000G))



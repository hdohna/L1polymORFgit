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
ScannedL1Ranges <- lapply(FileNames, function(x) scanBam(x))

# Count the number of reads mapped
NrMapped2L1 <- sapply(ScannedL1Ranges, function(x){
  sum(!is.na(x[[1]]$pos))
})
FilesWithReads <- FileNames[NrMapped2L1 > 0]

# Get aligned reads per peak
# R1 <- GRanges(seqnames = "L1HS_L1_Homo_sapiens", 
#               ranges = IRanges(start = 1, end = 6000))
# ReadsPerL1 <- lapply(FilesWithReads, function(x) {
#   Reads <- extractReads(x, R1)
# })
# 
# # Calculate a coverage matrix
# CoverMat <- t(sapply(ReadsPerL1, function(x){
#   Cov <- coverage(x)
#   as.vector(Cov[[1]])
# }))
CoverMat <- t(sapply(FilesWithReads, function(x) {
  RList <- scanBam(x)[[1]]
  CoverageFromReadList(RList, End = 6064)
}))
dim(CoverMat)

# Plot mean coverage
plot(colMeans(CoverMat), ylab = "Mean coverage", xlab = "Position on L1", type = "s")

# Get means and 95% quantiles 
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

# Full-length insertions
idx5P   <- which(CoverMat[,10] > 0)
Max3P   <- sapply(idx5P, function(x) max(which(CoverMat[x,] > 0)))
idxFull <- idx5P[Max3P >= 500]
x <- FilesWithReads[1]
plot(CoverMat[idxFull[1],])
FullL1Info <- t(sapply(FilesWithReads[idxFull], function(x){
  FPathSplit <- strsplit(x, "/")[[1]]
  FName      <- FPathSplit[length(FPathSplit)]
  FName      <- substr(FName, 1, nchar(FName) - 4)
  strsplit(FName, "_")[[1]]
}))
FullL1Info <- base::as.data.frame(FullL1Info, stringsAsFactors = F)
colnames(FullL1Info) <- c("chromosome", "idx")
FullL1Info$idx <- as.numeric(FullL1Info$idx)
IslGRanges_reduced[FullL1Info$idx]
overlapsAny(IslGRanges_reduced[FullL1Info$idx], L1CatalogGR)



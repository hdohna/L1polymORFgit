##############################################
#
# General description:
#
#   The following script reads a bam for 12878 and explores peak
#   overlap 

# Input:
#
#     BamFile: path to file that contains mapped reads
#     L1TableFileName: path to file that contains L1HS ranges in a table

# Output:
#   
#    : ...

##############################################

########################################
#                                      #
#  Source packages and set parameters  #
#                                      #
########################################

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(raster)
library(rgdal)
library(Matrix)
library(ShortRead)
library(csaw)
library(chipseq)
library(BSgenome.Hsapiens.UCSC.hg19)
length(BSgenome.Hsapiens.UCSC.hg19[['chr1']])

# Get all ranges of reads for per chromosome
Chromosomes <- paste("chr", 1:22, sep = "")
IslRange    <- 2*10^4
CalcNew_IslandPerCoverage <- T

# Files and folders
BamFile         <- "D:/L1polymORF/Data/NA12878-L1HS_S1_L001.dedup.unique.sorted.bam"
L1TableFileName <- "D:/L1polymORF/Data/L1_repeat_table.csv"

# Define a function to plot several coverage examples
PlotCoverExamples <- function(NrToPlot = 5, idxHighest = 0){
  Cols      <- rainbow(NrToPlot)
  plot(CoverMat[higestCov[1 + idxHighest], ], type = "l", ylab = 'Coverage', 
       xlab = "Genomic position", col = Cols[1])
  for (i in 2:NrToPlot){
    lines(CoverMat[higestCov[i + idxHighest], ], col = Cols[i])
  }
  segments(x0 = c(7000, 13000), y0 = c(0, 0), x1 = c(7000, 13000), 
           y1 = c(2000, 2000), col = "black")
}


#######################################
#                                     #
#    Read in data                     #
#                                     #
#######################################

# Read in table with L1 ranges
L1Table <- read.csv(L1TableFileName)

# Create names of subfamilies
repNChar    <- as.character(L1Table$repName)
SubFamiliesLookUp <- sapply(unique(repNChar), function(x){
  c(Name = x,
    SubFam = paste(strsplit(x, '[1-9]')[[1]][1:2], collapse = "1"))
})
NameMatch <- match(repNChar, SubFamiliesLookUp["Name",])
SubFamilies <- SubFamiliesLookUp["SubFam", NameMatch]

# Make some corrections to create a proper GRanges object with L1 Seqences
L1IRanges <- IRanges(start = L1Table$genoStart,
                     end = L1Table$genoEnd)
L1GRanges <- GRanges(seqnames = L1Table$genoName, ranges = L1IRanges,
                     strand = L1Table$strand)
GRanges_L1HSFL <- L1GRanges[width(L1GRanges) >= 6000 & SubFamilies == "L1HS"]
GRanges_L1PAFL <- L1GRanges[width(L1GRanges) >= 5900 & SubFamilies == "L1PA"]

#################################################
#                                               #
#    Determine coverage per L1HS range          #
#                                               #
#################################################

# Determine which ranges to analyze
CurrentRanges <- GRanges_L1HSFL

## Get ranges around each L1Hs and load intersecting reads
LargeL1Ranges <- resize(CurrentRanges, IslRange, fix = "center")
cat("*******   Getting reads overlapping with L1 ...   *******\n")
CoverMat <- matrix(0, nrow = length(CurrentRanges), ncol = IslRange)

# Fill in coverage
for (i in 1:length(LargeL1Ranges)){
  cat("Calculating coverage range", i, "of", length(LargeL1Ranges), "\n")
  Reads <- extractReads(bam.file = BamFile , LargeL1Ranges[i])
  Cov     <- GenomicRanges::coverage(Reads)
  CovVect <- as.vector(Cov[[1]])[start(LargeL1Ranges[i]):end(LargeL1Ranges[i])]
  CoverMat[i,] <- CovVect
}

# Re-arrange coverage according to strand
idxNegative <- which(strand(LargeL1Ranges) == '-')
for (i in idxNegative){
  CoverMat[i,] <- CoverMat[i, IslRange:1]
}

# Get means and 95% quantiles 
QuantileMat <- apply(CoverMat, 2, FUN = function(x) quantile(x, c(0.05, 0.5, 0.95)))
idxFw <- 1:ncol(CoverMat)
idxRv <- CoverMat:1
plot(QuantileMat[2,], type = "n", ylim = c(0, 150), ylab = 'Coverage', xlab = "Genomic position")
polygon(c(idxFw, idxRv), c(QuantileMat[1, idxFw], QuantileMat[3, idxRv]),
        col = "grey", border = NA)
lines(QuantileMat[2,], lwd = 1.2)
segments(x0 = c(7000, 13000), y0 = c(0, 0), x1 = c(7000, 13000), 
         y1 = c(200, 200), col = "red")
CreateDisplayPdf('D:/L1polymORF/Figures/L1HSCoverageUnique.pdf')

# Get the 5 L1s with the higest total coverage
TotalCov  <- rowSums(CoverMat)
higestCov <- order(TotalCov, decreasing = T)
PlotCoverExamples(idxHighest = 0)
CreateDisplayPdf('D:/L1polymORF/Figures/L1HSCoverageUnique_exampleHigh.pdf')
PlotCoverExamples(idxHighest = 200)
CreateDisplayPdf('D:/L1polymORF/Figures/L1HSCoverageUnique_exampleMedium.pdf')
PlotCoverExamples(idxHighest = 300)
CreateDisplayPdf('D:/L1polymORF/Figures/LHS1CoverageUnique_exampleLow.pdf')
hist(TotalCov, breaks = seq(0, 4*10^5, 10000))
TotalCov[higestCov[200]]

# Overlay all coverage plot from coverage rank 100-200
PlotCoverExamples(NrToPlot = 100, idxHighest = 100)

# Get reads for L1 ranges
param <- ScanBamParam(which = LargeL1Ranges, what = scanBamWhat())
ScannedReads <- scanBam(file = BamFile, 
                         param = param)

# Create a matrix that counts paired reads
i <- 1
PairMatrix <- matrix(0, 200, 200)
for (i in 1:length(ScannedReads)){
  ReadData <- ScannedReads[[i]]
  if (as.vector(strand(LargeL1Ranges)[i]) == '+'){
    relPos  <- ReadData$pos  - start(LargeL1Ranges)[i]
    relMPos <- ReadData$mpos - start(LargeL1Ranges)[i]
  } else {
    relPos  <- end(LargeL1Ranges)[i] - ReadData$pos
    relMPos <- end(LargeL1Ranges)[i] - ReadData$mpos
  }
  relMPos[relMPos < 1 | relMPos > IslRange] <- NA
  PosInt  <- cut(relPos,  breaks = seq(0, IslRange, 100))
  MPosInt <- cut(relMPos, breaks = seq(0, IslRange, 100))
  CountTable <- table(PosInt, MPosInt)
  PairMatrix <- PairMatrix + CountTable
}
dim(CountTable)
graphics::image(PairMatrix, col = topo.colors(10))
L1Borders <- c(7000, 13000)/20000
segments(x0 = c(0, 0, L1Borders), y0 = c(L1Borders, 0, 0), x1 = c(200, 200, L1Borders), 
         y1 = c(L1Borders, 200, 200), col = "red")
lines(1:ncol(QuantileMat) / 20000, QuantileMat[2,] / (2*max(QuantileMat[2,])), col = "red")
CreateDisplayPdf('D:/L1polymORF/Figures/L1HSPairPlot.pdf')

#################################################
#                                               #
#    Determine coverage per L1PA range          #
#                                               #
#################################################

# Determine which ranges to analyze
CurrentRanges <- GRanges_L1PAFL

## Get ranges around each L1Hs and load intersecting reads
LargeL1Ranges <- resize(CurrentRanges, IslRange, fix = "center")
cat("*******   Getting reads overlapping with L1 ...   *******\n")
CoverMat <- matrix(0, nrow = length(CurrentRanges), ncol = IslRange)

# Fill in coverage
for (i in 1:length(LargeL1Ranges)){
  cat("Calculating coverage range", i, "of", length(LargeL1Ranges), "\n")
  Reads   <- extractReads(bam.file = BamFile , LargeL1Ranges[i])
  Cov     <- GenomicRanges::coverage(Reads)
  CovVect <- as.vector(Cov[[1]])[start(LargeL1Ranges[i]):end(LargeL1Ranges[i])]
  CoverMat[i,] <- CovVect
}

# Re-arrange coverage according to strand
idxNegative <- which(strand(LargeL1Ranges) == '-')
for (i in idxNegative){
  CoverMat[i,] <- CoverMat[i, IslRange:1]
}

# Get means and 95% quantiles 
QuantileMat <- apply(CoverMat, 2, FUN = function(x) quantile(x, c(0.05, 0.5, 0.95)))
idxFw <- 1:ncol(CoverMat)
idxRv <- ncol(CoverMat):1
plot(QuantileMat[2,], type = "n", ylim = c(0, 100), ylab = 'Coverage', xlab = "Genomic position")
polygon(c(idxFw, idxRv), c(QuantileMat[1, idxFw], QuantileMat[3, idxRv]),
        col = "grey", border = NA)
lines(QuantileMat[2,], lwd = 1.2)
segments(x0 = c(7000, 13000), y0 = c(0, 0), x1 = c(7000, 13000), 
         y1 = c(2000, 2000), col = "red")
CreateDisplayPdf('D:/L1polymORF/Figures/L1PACoverageUnique.pdf')

# Get the 5 L1s with the higest total coverage
TotalCov  <- rowSums(CoverMat)
higestCov <- order(TotalCov, decreasing = T)
PlotCoverExamples(idxHighest = 0)
CreateDisplayPdf('D:/L1polymORF/Figures/L1PACoverageUnique_exampleHigh.pdf')
PlotCoverExamples(idxHighest = 3000)
CreateDisplayPdf('D:/L1polymORF/Figures/L1PACoverageUnique_exampleMedium.pdf')
PlotCoverExamples(idxHighest = 6500)
CreateDisplayPdf('D:/L1polymORF/Figures/LHPACoverageUnique_exampleLow.pdf')
hist(TotalCov, breaks = seq(0, 4*10^5, 10000))
TotalCov[higestCov[200]]
nrow(CoverMat)

# Get reads for L1 ranges
param <- ScanBamParam(which = LargeL1Ranges, what = scanBamWhat())
ScannedReads <- scanBam(file = BamFile, 
                        param = param)

# Create a matrix that counts paired reads
PairMatrix <- matrix(0, 200, 200)
for (i in 1:length(ScannedReads)){
  ReadData <- ScannedReads[[i]]
  if (as.vector(strand(LargeL1Ranges)[i]) == '+'){
    relPos  <- ReadData$pos  - start(LargeL1Ranges)[i]
    relMPos <- ReadData$mpos - start(LargeL1Ranges)[i]
  } else {
    relPos  <- end(LargeL1Ranges)[i] - ReadData$pos
    relMPos <- end(LargeL1Ranges)[i] - ReadData$mpos
  }
  relMPos[relMPos < 1 | relMPos > IslRange] <- NA
  PosInt  <- cut(relPos,  breaks = seq(0, IslRange, 100))
  MPosInt <- cut(relMPos, breaks = seq(0, IslRange, 100))
  CountTable <- table(PosInt, MPosInt)
  PairMatrix <- PairMatrix + CountTable
}
dim(CountTable)
graphics::image(PairMatrix, col = topo.colors(10))
L1Borders <- c(7000, 13000)/20000
segments(x0 = c(0, 0, L1Borders), y0 = c(L1Borders, 0, 0), x1 = c(200, 200, L1Borders), 
         y1 = c(L1Borders, 200, 200), col = "red")
lines(1:ncol(QuantileMat) / 20000, QuantileMat[2,] / (2*max(QuantileMat[2,])), col = "red")
CreateDisplayPdf('D:/L1polymORF/Figures/L1PAPairPlot.pdf')


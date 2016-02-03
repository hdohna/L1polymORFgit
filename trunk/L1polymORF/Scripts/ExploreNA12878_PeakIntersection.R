##############################################
#
# General description:
#
#   The following script reads a bam for 12878 and explores peak
#   overlap 

# Input:
#
#     BamFile: path to file that contains mapped reads
#     L1HSTableFileName: path to file that contains L1HS ranges in a table

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
BamFile            <- "D:/L1polymORF/Data/NA12878-L1HS_S1_L001.dedup.unique.sorted.bam"
L1HSTableFileName  <- "D:/L1polymORF/Data/L1HS_repeat_table.csv"
scanBamHeader(BamFile)

#######################################
#                                     #
#    Read in data                     #
#                                     #
#######################################

# Read in table with L1 ranges
L1HSTable <- read.csv(L1HSTableFileName)

# Make some corrections to create a proper GRanges object with L1 Seqences
L1IRanges <- IRanges(start = L1HSTable$genoStart,
                     end = L1HSTable$genoEnd)
L1GRanges <- GRanges(seqnames = L1HSTable$genoName, ranges = L1IRanges,
                     strand = L1HSTable$strand)
L1GRanges <- L1GRanges[width(L1GRanges) >= 6000]


# Get ranges around each L1Hs and load intersecting reads
LargeL1Ranges <- resize(L1GRanges, IslRange, fix = "center")
cat("*******   Getting reads overlapping with L1 ...   *******\n")
CoverMat <- matrix(0, nrow = length(L1GRanges), ncol = IslRange)

# Fill in coverage
for (i in 1:length(LargeL1Ranges)){
  print(i)
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
idxFw <- 1:length(CovMean)
idxRv <- length(CovMean):1
plot(QuantileMat[2,], type = "n", ylim = c(0, 150), ylab = 'Coverage', xlab = "Genomic position")
polygon(c(idxFw, idxRv), c(QuantileMat[1, idxFw], QuantileMat[3, idxRv]),
        col = "grey", border = NA)
lines(QuantileMat[2,], lwd = 1.2)
segments(x0 = c(7000, 13000), y0 = c(0, 0), x1 = c(7000, 13000), 
         y1 = c(200, 200), col = "red")
CreateDisplayPdf('D:/L1polymORF/Figures/L1CoverageUnique.pdf')

# Get the 5 L1s with the higest total coverage
TotalCov  <- rowSums(CoverMat)
higestCov <- order(TotalCov, decreasing = T)
PlotCoverExamples <- function(NrToPlot = 5, idxHighest = 0){
  Cols      <- rainbow(NrToPlot)
  plot(CoverMat[higestCov[1 + idxHighest], ], type = "l", ylab = 'Coverage', 
       xlab = "Genomic position", col = Cols[1])
  for (i in 2:NrToPlot){
    lines(CoverMat[higestCov[i + idxHighest], ], col = Cols[i])
  }
  segments(x0 = c(7000, 13000), y0 = c(0, 0), x1 = c(7000, 13000), 
           y1 = c(200, 200), col = "black")
}
PlotCoverExamples(idxHighest = 0)
CreateDisplayPdf('D:/L1polymORF/Figures/L1CoverageUnique_exampleHigh.pdf')
PlotCoverExamples(idxHighest = 200)
CreateDisplayPdf('D:/L1polymORF/Figures/L1CoverageUnique_exampleMedium.pdf')
PlotCoverExamples(idxHighest = 300)
CreateDisplayPdf('D:/L1polymORF/Figures/L1CoverageUnique_exampleLow.pdf')
hist(TotalCov, breaks = seq(0, 4*10^5, 10000))
TotalCov[higestCov[200]]

PlotCoverExamples(NrToPlot = 100, idxHighest = 100)

param <- ScanBamParam(which = LargeL1Ranges, what = scanBamWhat())
ScannedReads <- scanBam(file = BamFile, 
                         param = param)

i <- 1
PairMatrix <- matrix(0, 200, 200)
for (i in 1:length(ScannedReads)){
  ReadData <- ScannedReads[[i]]
  relPos  <- ReadData$pos  - start(LargeL1Ranges)[i]
  relMPos <- ReadData$mpos - start(LargeL1Ranges)[i]
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
CreateDisplayPdf('D:/L1polymORF/Figures/PairPlot.pdf')


##############################################
#
# General description:
#
#   The following script reads a bam for 12878 and explores peak
#   overlap and reads in pacbio data for the same ranges

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
BamFile            <- "D:/L1polymORF/Data/NA12878-L1HS_S1_L001.dedup.unique.sorted.bam"
L1HSTableFileName  <- "D:/L1polymORF/Data/L1HS_repeat_table.csv"
L1HSTableFileName19  <- "D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv"
AlignFileName      <- "D:/L1polymORF/Data/L1Sequences_reptab_alignedMAFFT"

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
L1Table <- read.csv(L1HSTableFileName)

# Create GRanges object with L1 Seqences
L1IRanges <- IRanges(start = L1Table$genoStart,
                     end = L1Table$genoEnd)
L1GRanges <- GRanges(seqnames = L1Table$genoName, ranges = L1IRanges,
                     strand = L1Table$strand)
L1GRanges <- L1GRanges[width(L1GRanges) >= 6000] 

#################################################
#                                               #
#    Determine coverage per L1HS range          #
#                                               #
#################################################

# Determine which ranges to analyze
CurrentRanges <- L1GRanges

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
idxRv <- ncol(CoverMat):1

#################################################
#                                               #
#    Load Pacbio reads from the same ranges     #
#                                               #
#################################################

# Load results
load("D:/L1polymORF/Data/ReadsPerL1.RData")
NrMapped <- sapply(ScannedReads, function(x) sum(!is.na(x$pos)))
hist(NrMapped, breaks= seq(0, 20000, 10),xlim = c(0, 1000))

# Deterimine the coverage within the L1 range
CoverInL1  <- rowSums(CoverMat[,7000:13000])
CoverOutL1 <- rowSums(CoverMat[,-c(7000:13000)])
CoverInL1[which(NrMapped == 0)]
CoverOutL1/CoverInL1
hist(CoverInL1, breaks = seq(0, 290000, 1000))
min(CoverInL1)

plot(NrMapped, CoverInL1, xlim = c(0, 1000))
plot(NrMapped, CoverOutL1/CoverInL1, xlim = c(0, 1000),
     ylim = c(0, 1))

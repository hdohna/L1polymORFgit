# The script below subsets creates GRanges with L1s from the 1000 Genome data
# Phase 3, available at ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
# The granges are created from a L1 table created by the script 
# 'Create_1000G_L1Table.R'

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(GenomicRanges)
library(rtracklayer)

# Specify file paths
GROutputPath <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
L1RefRangePath <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
load(L1RefRangePath)

# Specify parameters
NrInfoCols <- 9

# Read in table with L1s from 1000 genome data
cat("Reading L1 table\n")
L1_1000G <- read.table("D:/L1polymORF/Data/L1_1000G_withGenoNum")
colnames(L1_1000G)[1:10]

# Get the sample columns
SampleColumns <- colnames(L1_1000G)[(NrInfoCols + 1):(ncol(L1_1000G) - 1)]
SampleColumns[length(SampleColumns)]

# Add additional information to L1_1000G
L1_1000G$chromosome <- paste("chr", L1_1000G$CHROM, sep = "")
L1_1000G$Frequency  <- rowSums(L1_1000G[,SampleColumns]) / 2 / length(SampleColumns)
L1InfoCols <- data.frame(t(sapply(1:nrow(L1_1000G), function(i){
  x <- as.character(L1_1000G$INFO[i])
  InfoSplit   <- strsplit(x, ";")[[1]]
  L1Info      <- grep("MEINFO=", InfoSplit, value = T)
  if(length(L1Info) > 0){
    L1InfoSplit <- strsplit(L1Info, ",")[[1]]
  } else {
    L1InfoSplit <- rep(NA, 4)
  }
  c(L1Start = L1InfoSplit[2], L1End = L1InfoSplit[3], L1Strand = L1InfoSplit[4])
})))
L1_1000G <- cbind(L1_1000G, L1InfoCols)
L1_1000G_reduced <- L1_1000G[,
     c("chromosome", "POS", "Frequency", "InsLength", "L1Start", "L1End", 
       "L1Strand")]

# Create genomic ranges of of L1 table
L1_1000G_GR_hg19 <- makeGRangesFromDataFrame(L1_1000G_reduced, 
                                       keep.extra.columns = T,
                                       start.field="POS",
                                       end.field ="POS")
# Create a subset for NA12878 
L1_1000G_GR_hg19_NA12878 <- L1_1000G_GR_hg19[L1_1000G$NA12878 > 0]

# Lift over to other genome builds
cat("Lift over to hg38\n")
L1_1000G_GRList_hg38 <- UniqueLiftover(L1_1000G_GR_hg19,
   ChainFilePath = "D:/L1polymORF/Data/hg19ToHg38.over.chain")

# Save genomic ranges
cat("Saving genomic ranges\n")
save(list = c("L1_1000G_GR_hg19", "L1_1000G_GR_hg19_NA12878", "L1_1000G_GRList_hg38", 
              "L1_1000G_reduced", "SampleColumns", "L1_1000G"), 
     file = GROutputPath)

# Get for each 1000 genome L1 the closest distance to a reference L1
ClostestDist <- Dist2Closest(L1_1000G_GR_hg19, L1GRanges)
hist(ClostestDist, breaks = seq(0, 60000, 250))

# Plot frequency histogram of full-length and 
hist(L1_1000G$Frequency[L1_1000G$InsLength > 6000], breaks = seq(0, 1, 0.01))
hist(L1_1000G$Frequency[L1_1000G$InsLength < 6000], breaks = seq(0, 1, 0.001))
mean(L1_1000G$Frequency[L1_1000G$InsLength > 6000])
mean(L1_1000G$Frequency[L1_1000G$InsLength < 6000])
sum(L1_1000G$InsLength > 6010)
sum(L1_1000G$Frequency == 0)
L1_1000G_GR_hg19@elementMetadata@listData$Frequency

# Get the sum of full-length L1 per individual
blnFull <- L1_1000G$InsLength > 6000
SumFull <- sapply(SampleColumns, function(x) sum(L1_1000G[blnFull, x]))
ExpFull <- 2 * sum(L1_1000G$Frequency[blnFull])
mean(SumFull < ExpFull)

# Get the sum of fragment L1 per individual
SumFragm <- sapply(SampleColumns, function(x) sum(L1_1000G[!blnFull, x]))
ExpFragm <- 2 * sum(L1_1000G$Frequency[!blnFull])
mean(SumFragm < ExpFragm)

##########
# Process L1 catalog
##########

cat("Comparing to L1 catalog\n")

# Read in table with known L1 
L1Catalogue <- read.csv("D:/L1polymORF/Data/L1CatalogExtended.csv", as.is = T)
L1Catalogue$Allele[is.na(L1Catalogue$Allele)] <- 1

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
blnInRef1         <- (L1Catalogue$end_HG38 - L1Catalogue$start_HG38) > 6000 
L1CatalogL1Mapped <- L1Catalogue[blnL1Mapped & blnAllele1,]

# Create genomic ranges for catalog L1
L1CatalogGR <- GRanges(seqnames = L1CatalogL1Mapped$Chromosome,
                       ranges = IRanges(start = pmin(L1CatalogL1Mapped$start_HG38,
                                                     L1CatalogL1Mapped$end_HG38),
                                        end = pmax(L1CatalogL1Mapped$start_HG38,
                                                   L1CatalogL1Mapped$end_HG38)),
                       strand = L1CatalogL1Mapped$strand_L1toRef)

DistCat2_1000G <- Dist2Closest(L1CatalogGR, L1_1000G_GRList_hg38$LiftedRanges)

# Get indices of 1000 Genome and catalog elements that match
idx1000G <- nearest(L1CatalogGR, L1_1000G_GRList_hg38$LiftedRanges)
L1CatalogMatch1000G <- L1CatalogL1Mapped[DistCat2_1000G < 100, ]
idx1000GMatchCat    <- idx1000G[DistCat2_1000G < 100]

# Get rows of L1_1000G table that can be matched to catalog
L1_1000G_match <- L1_1000G[idx1000GMatchCat, ]
ExpectedAct <- 2 * L1CatalogMatch1000G$ActivityNum %*% L1_1000G_match$Frequency
ObservedAct <- sapply(SampleColumns, function(x){
  L1CatalogMatch1000G$ActivityNum %*% L1_1000G_match[,x]})
hist(ObservedAct)
mean(ObservedAct)

# Plot activity vs frequency
plot(L1CatalogMatch1000G$ActivityNum, L1_1000G_match$Frequency)
cor(L1CatalogMatch1000G$ActivityNum, L1_1000G_match$Frequency)
cor.test(L1CatalogMatch1000G$ActivityNum, L1_1000G_match$Frequency)

# Sample individual activity sums
NrSamples <- 100000
SampledActSums <- sapply(1:NrSamples, function(x){
  rbinom(nrow(L1_1000G_match), 2, L1_1000G_match$Frequency) %*% 
    L1CatalogMatch1000G$ActivityNum
})
# cat("Sampling quantiles of activity sums\n")
# SamplQuantList <- SampleQuantiles(SampledActSums, length(SampleColumns), 
#                                   NrSamples = 10000, QuantV = seq(0, 1, 0.0001),
#                                   LowerQ = 0.005, UpperQ = 0.995)
# QSMat <- SamplQuantList$QMat
# idxF <- 1:ncol(QSMat)
# idxR <- ncol(QSMat):1
# mean(SampledActSums)
# QQ1  <- qqplot(SampledActSums, ObservedAct, plot.it = F)
# plot(QQ1$x, QQ1$y, xlab = "Sampled activity sums", 
#      ylab = "Observed activity sums")
# MedianMatch1 <- sapply(QQ1$x, function(z) which.min(abs(z - QSMat[1, ])))
# MedianMatch2 <- sapply(QQ1$x, function(z) which.min(abs(z + 1 - QSMat[1, ])))
# MedianMatch3 <- sapply(QQ1$x, function(z) which.min(abs(z - 1 - QSMat[1, ])))
# polygon(QSMat[1, c(idxF, idxR)], c(QSMat[2, ], QSMat[3, idxR]), 
#         col = "grey", border = NA)
# points(QQ1$x, QQ1$y)
# lines(c(0, 1000), c(0, 1000))
# blnOutside <- (QQ1$y < QSMat[2, MedianMatch1] | QQ1$y > QSMat[3, MedianMatch1]) &
#   (QQ1$y < QSMat[2, MedianMatch2] | QQ1$y > QSMat[3, MedianMatch2]) &
#   (QQ1$y < QSMat[2, MedianMatch3] | QQ1$y > QSMat[3, MedianMatch3])
# sum(blnOutside)
# points(QQ1$x[blnOutside],  QQ1$y[blnOutside], col = "red")
#   
# CreateDisplayPdf('D:/L1polymORF/Figures/L1ActivitySums.pdf', 
#                  PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

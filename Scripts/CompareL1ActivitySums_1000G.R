# The script below reads GRanges with L1s from the 1000 Genome data
# Phase 3, available at ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
# The script takes matches L1 insertions to 
# The Granges were created by the script 'Create_1000G_L1GRanges.R'

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(GenomicRanges)
library(rtracklayer)

# Specify file paths
L1GRPath       <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
L1RefRangePath <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'

# Specify parameters
NrInfoCols <- 9

# Load previously generated objects
load(L1RefRangePath)
load(L1GRPath)


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
length(ObservedAct)

# Plot activity vs frequency
plot(L1CatalogMatch1000G$ActivityNum, L1_1000G_match$Frequency)
cor(L1CatalogMatch1000G$ActivityNum, L1_1000G_match$Frequency)
cor.test(L1CatalogMatch1000G$ActivityNum, L1_1000G_match$Frequency)

# Sample individual activity sums
NrSamples <- 1000000
SampledActSums <- sapply(1:NrSamples, function(x){
  rbinom(nrow(L1_1000G_match), 2, L1_1000G_match$Frequency) %*% 
    L1CatalogMatch1000G$ActivityNum
})

# Get samples of variances
StartVals <- seq(1, NrSamples, length(ObservedAct))
SampledVars <- sapply(StartVals[-length(StartVals)], function(x){
  var(SampledActSums[x:(x+length(ObservedAct))])
})
par(mar = c(5, 4, 4, 2) + 0.1)
hist(SampledVars, main = "", xlab = "Variance in activity sum")
segments(var(ObservedAct), y0 = 0, y1 = 50, col = "red")
PVar <- mean(SampledVars <= var(ObservedAct))
CreateDisplayPdf('D:/L1polymORF/Figures/L1VarianceOfActivitySums.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')


cat("Sampling quantiles of activity sums\n")
SamplQuantList <- SampleQuantiles(SampledActSums, length(SampleColumns), 
                                  NrSamples = 10000, QuantV = seq(0, 1, 0.0001),
                                  LowerQ = 0.005, UpperQ = 0.995)
QSMat <- SamplQuantList$QMat
idxF <- 1:ncol(QSMat)
idxR <- ncol(QSMat):1
mean(SampledActSums)
QQ1  <- qqplot(SampledActSums, ObservedAct, plot.it = F)
plot(QQ1$x, QQ1$y, xlab = "Sampled activity sums", 
     ylab = "Observed activity sums")
MedianMatch1 <- sapply(QQ1$x, function(z) which.min(abs(z - QSMat[1, ])))
MedianMatch2 <- sapply(QQ1$x, function(z) which.min(abs(z + 1 - QSMat[1, ])))
MedianMatch3 <- sapply(QQ1$x, function(z) which.min(abs(z - 1 - QSMat[1, ])))
polygon(QSMat[1, c(idxF, idxR)], c(QSMat[2, ], QSMat[3, idxR]), 
        col = "grey", border = NA)
points(QQ1$x, QQ1$y)
lines(c(0, 1000), c(0, 1000))
# blnOutside <- (QQ1$y < QSMat[2, MedianMatch1] | QQ1$y > QSMat[3, MedianMatch1]) &
#   (QQ1$y < QSMat[2, MedianMatch2] | QQ1$y > QSMat[3, MedianMatch2]) &
#   (QQ1$y < QSMat[2, MedianMatch3] | QQ1$y > QSMat[3, MedianMatch3])
blnOutside <- (QQ1$y < QSMat[2, MedianMatch1] | QQ1$y > QSMat[3, MedianMatch1])
sum(blnOutside)
points(QQ1$x[blnOutside],  QQ1$y[blnOutside], col = "red")
  
CreateDisplayPdf('D:/L1polymORF/Figures/L1ActivitySums.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

# Compare histogram of observed activity sum with histograms of simulated sums
BinWidth <- 25
HistBreaks <- seq(0, 2000, BinWidth)
Obs = ObservedAct
ObsDens <- hist(ObservedAct, HistBreaks, plot = F)$density
SampledDens <- sapply(StartVals[-length(StartVals)], function(x){
  hist(SampledActSums[x:(x+length(ObservedAct))], HistBreaks, plot = F)$density
})
dim(SampledDens)
SampledDensMean <- rowMeans(SampledDens)
SampledDensQMat <- apply(SampledDens, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))

plot(c(0, HistBreaks[-length(HistBreaks)]), c(0, ObsDens), 
     type = "s", xlim = c(0, 600))
lines(c(0, HistBreaks[-length(HistBreaks)]), c(0, SampledDensMean),
      type = "s", col = "red")

plot(HistBreaks[-length(HistBreaks)], (ObsDens - SampledDensMean) / SampledDensMean,
     xlim = c(0, 1000))
lines(c(0, 5000), c(0, 0), col = "red", lty = 2)
mean(ObservedAct)

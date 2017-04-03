# The following script reads in a Hi-C data and determines how much
# different L1 insertions interact with other genomic regions

library(GenomicRanges)
library(rtracklayer)
library(Matrix)
library(irlba)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.r')

############
#  Load results
############

load("D:/L1polymORF/Data/HiCData/10kb_L1summary.Rdata")

# Plot 
boxplot(NormEOsum ~ L1type, HicByL1Type)

# Get HiC values for full-length and fragment L1
blnFragm    <- HicByL1Type$L1type == "fragm"
HiCL1fragm  <- HicByL1Type$NormEOsum[blnFragm]
MeanFragm    <- mean(HiCL1fragm)

blnFull     <- HicByL1Type$L1type == "full"
NrFull      <- sum(blnFull)
HiCL1Full   <- HicByL1Type$NormEOsum[blnFull]
MeanFull    <- mean(HiCL1Full)

# Sample HiC values from fragments
SQuant     <- SampleQuantiles(HiCL1fragm, SubsetSize = NrFull, 
                              NrSamples = 10^4)
  
# Plot histogram of means sampled from fragments against the oberved mean of 
# full-length L1
hist(SQuant$SampleMeans, xlab = "Normalized EO")
segments(x0 = MeanFull, y0 = 0, x1 = MeanFull, y1 = 1000, lwd = 3)

# Get p-value from samples and t-test
mean(SQuant$SampleMeans < MeanFull)
t.test(NormEOsum ~ L1type, HicByL1Type)

# P-value per chromosoeme
pValPerChrom

qqplot(HiCL1fragm, HicByL1Type$NormEOsum[blnFull])
lines(c(0, 10^10), c(0, 10^10))
max(HiCL1fragm)

# Plot histogram of HiC values per fragment and full-length L1
HistFragm <- hist(HiCL1fragm, breaks = seq(0, 85000, 5000), plot = F)
HistFull  <- hist(HiCL1Full, breaks = seq(0, 85000, 5000), plot = F)
plot(HistFragm$mids, HistFragm$density, type = "s", col = "red", ylim = c(0, 0.0001))
lines(HistFull$mids, HistFull$density, type = "s", col = "blue")
segments(x0 = MeanFull, y0 = 0, x1 = MeanFull, y1 = 0.0002, lwd = 3, col = "blue")
segments(x0 = MeanFragm, y0 = 0, x1 = MeanFragm, y1 = 0.0002, lwd = 3, col = "red")

# Plot number of fragment and full-length L1 per chromosome
L1PerChrom <- table(HicByL1Type$L1type, HicByL1Type$chromosome)
ChrMatch   <- match(colnames(L1PerChrom), names(ChromLengthsHg19))
plot(ChromLengthsHg19[ChrMatch], L1PerChrom["fragm", ])
text(ChromLengthsHg19[ChrMatch], L1PerChrom["fragm", ], names(ChromLengthsHg19[ChrMatch]))

barplot(L1PerChrom["fragm", ] / ChromLengthsHg19[ChrMatch], las = 2)
barplot(L1PerChrom["full", ]/ChromLengthsHg19[ChrMatch], las = 2)

# Determine correlation between L1 fragment length and HiC value
blnFragm <- HicByL1Type$L1type == "fragm"
L1Width  <- HicByL1Type$end - HicByL1Type$start
cor(L1Width[blnFragm], HicByL1Type$NormEOsum[blnFragm])
cor.test(L1Width[blnFragm], HicByL1Type$NormEOsum[blnFragm])
plot(L1Width[blnFragm], HicByL1Type$NormEOsum[blnFragm])
HiCVsLengthSmu <- supsmu(L1Width[blnFragm], HicByL1Type$NormEOsum[blnFragm])
lines(HiCVsLengthSmu$x, HiCVsLengthSmu$y, col = "red")

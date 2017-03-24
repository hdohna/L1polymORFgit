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

load("D:/L1polymORF/Data/HiCData/50kb_L1summary.Rdata")
HicByL1Type$L1type
boxplot(NormEOsum ~ L1type, HicByL1Type)
t.test(NormEOsum ~ L1type, HicByL1Type)
pValPerChrom


# Sample HiC values
HiCL1fragm <- HicByL1Type$NormEOsum[HicByL1Type$L1type == "fragm"]
blnFull    <- HicByL1Type$L1type == "full"
NrFull     <- sum(blnFull)
MeanFull   <- mean(HicByL1Type$NormEOsum[blnFull])
SQuant     <- SampleQuantiles(HiCL1fragm, SubsetSize = NrFull)
  
# Plot histogram of means sampled from fragments against the oberved mean of 
# full-length L1
hist(SQuant$SampleMeans, xlab = "Normalized EO")
segments(x0 = MeanFull, y0 = 0, x1 = MeanFull, y1 = 300, lwd = 3)
mean(SQuant$SampleMeans < MeanFull)

qqplot(HiCL1fragm, HicByL1Type$NormEOsum[blnFull])
lines(c(0, 10^10), c(0, 10^10))

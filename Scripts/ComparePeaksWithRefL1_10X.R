# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(ape)
library(ShortRead)
library(rtracklayer)
library(Rsamtools)
library(csaw)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

ComparePeaksWithRefL1(
  BamFile = 'D:/L1polymORF/Data/NA12878_capt10X_aln2hg19_NonRefL1filtered.bam',
  OutBamFileFullLengthL1 = NULL,
  L1Ranges = "D:/L1polymORF/Data/L1RefRanges_hg19.Rdata",
  MinMaxCover = 10,    # minimum maximum coverage to be called a peak 
  MinGap      = 1000,
  NrChromPieces = 1,
  MinDist2L1  = 3*10^4, # minimum distance to L1 to be called a peak 
  OutFile = "D:/L1polymORF/Data/Analyzed10X1Ranges_hg19.RData",
  EndList = NULL,
  blnFilterOverlap = F
)
load("D:/L1polymORF/Data/Analyzed10X1Ranges_hg19.RData")
hist(maxCover, xlim = c(0, 100), breaks = 0:300000)
sum(maxCover == 1)
hist(width(IslGRanges_reduced), breaks = seq(0, 200000, 1000), xlim = c(0, 50000))
plot(width(IslGRanges_reduced), maxCover, ylim = c(0, 1000))
cor(width(IslGRanges_reduced), maxCover)

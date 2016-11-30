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
  BamFile = 'D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_subreads_hg19withL1.sorted.bam',
  OutBamFileFullLengthL1 = NULL,
  L1Ranges = "D:/L1polymORF/Data/L1RefRanges_hg19.Rdata",
  MinMaxCover = 2,    # minimum maximum coverage to be called a peak 
  MinGap      = 1000,
  NrChromPieces = 1,
  MinDist2L1  = 3*10^4, # minimum distance to L1 to be called a peak 
  OutFile = "D:/L1polymORF/Data/AnalyzedPacBioL1Ranges_MinGap1000.RData",
  EndList = NULL,
  blnFilterOverlap = F
)
  
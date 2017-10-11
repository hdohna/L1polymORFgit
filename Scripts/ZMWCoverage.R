# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(seqinr)
library(ShortRead)
library(Rsamtools)
library(rtracklayer)
library(csaw)
library(GenomicRanges)
#library(BSgenome.Hsapiens.UCSC.hg19)

# Create a read list object
RL    <- scanBam("D:/L1polymORF/Data/BZ_NonRef/chr16_4873.bam")[[1]]
NotNa <- !is.na(RL$pos)
RL    <- lapply(RL, function(x) x[NotNa])
ChrLength <- 6064

ZMWCover <- ZMWCoverage(RL, ChrLength = 6064)
slice(ZMWCover, lower = 5)


# Create a read list object
ScanRange <- GRanges("L1HS_L1_Homo_sapiens", IRanges(start = 3816, end = 3816))
ScanPars <- ScanBamParam(what = scanBamWhat(), which = ScanRange)
RLTest    <- scanBam("D:/L1polymORF/Data/BZ_NonRef/chr16_4873.bam", param = ScanPars)[[1]]
NotNa     <- !is.na(RLTest$pos)
RLTest  <- lapply(RLTest, function(x) x[NotNa])
extractReads("D:/L1polymORF/Data/BZ_NonRef/chr16_4873.bam", ScanRange, as.reads = T)

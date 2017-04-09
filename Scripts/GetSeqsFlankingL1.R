# The script below reads data on L1 ranges and saves flanking sequences

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(smooth)
library(ShortRead)
library(Rsamtools)
library(rtracklayer)
library(csaw)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

# Load data on L1 coverage
load('D:/L1polymORF/Data/L1_NA12878_PacBio_Coverage.RData')

# Read table with L1 insertion info (created in script 'AnalyzeReadsMapped2L1.R')
FullL1Info <- read.csv("D:/L1polymORF/Data/L1InsertionInfo.csv", row.names = 1)

# Subset to get full-length insertions
idxFull    <- which(sapply(1:nrow(CoverMat), function(x) 
  all(CoverMat[x, 50:6040] > 0)))
NamesFull  <- rownames(CoverMat)[idxFull]
blnFullIns <- rownames(FullL1Info) %in% NamesFull
FullL1InfoSubset <- FullL1Info[blnFullIns, ]

# Get genomic ranges and and sequences of flanks
L1LeftFlanks <- GRanges(seqnames = FullL1InfoSubset$chromosome, 
   ranges = IRanges(end = FullL1InfoSubset$L1InsertionPosition.median,
                                         width = 1000))
L1RightFlanks <- GRanges(seqnames = FullL1InfoSubset$chromosome, 
   ranges = IRanges(start = FullL1InfoSubset$L1InsertionPosition.median,
                                         width = 1000))
L1LeftFlankSeq  <- getSeq(BSgenome.Hsapiens.UCSC.hg19, L1LeftFlanks)
L1RightFlankSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, L1RightFlanks)
names(L1LeftFlankSeq) <- paste(seqnames(L1LeftFlanks), end(L1LeftFlanks), 
                               "LeftFlank", sep = "_")
names(L1RightFlankSeq) <- paste(seqnames(L1RightFlanks), 
   start(L1RightFlanks),  "RightFlank", sep = "_")

# Write out flanking sequences as fasta file
writeFasta(L1LeftFlankSeq,  "D:/L1polymORF/Data/L1LeftFlank.fasta")
writeFasta(L1RightFlankSeq, "D:/L1polymORF/Data/L1RightFlank.fasta")

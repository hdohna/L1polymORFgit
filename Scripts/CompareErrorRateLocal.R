# This script summarizes the number of mismatches with reference sequences 
# per phred score

# Load packages
library(seqinr)
library(ShortRead)
library(Rsamtools)
library(rtracklayer)
library(csaw)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

# Paths (local computer)
StartPath <- 'D:/L1polymORFgit/Scripts/_Start_L1polymORF.R'
DataPath  <- "D:/L1polymORF/Data/ErrorComparer.RData"

# Source start script
source(StartPath)

# Load data
load(DataPath)

# Calculate means
ErrorHiFi
ErrorHiFiBWA
ErrorNormal
ErrorNormalBWA
ErrorRefHiFi      
ErrorRefHiFiBWA   
ErrorRefNormal    
ErrorRefNormalBWA 
mean(ErrorHiFi,         na.rm = T)
mean(ErrorHiFiBWA,      na.rm = T)
mean(ErrorNormal,       na.rm = T)
mean(ErrorNormalBWA,    na.rm = T)
mean(ErrorRefHiFi,      na.rm = T)
mean(ErrorRefHiFiBWA,   na.rm = T)
mean(ErrorRefNormal,    na.rm = T)
mean(ErrorRefNormalBWA, na.rm = T)

sum(!is.na(ErrorHiFi))
sum(!is.na(ErrorHiFiBWA))
sum(!is.na(ErrorNormal))
sum(!is.na(ErrorNormalBWA))
sum(!is.na(ErrorHiFi))
sum(!is.na(ErrorHiFi))

# Find overlaps between the two genomic ranges and 
GRHiFi      <- unlist(GRangesList(GRHiFi))
GRNormalBWA <- unlist(GRangesList(GRNormalBWA))
ROverlap    <- findOverlaps(GRHiFi, GRNormalBWA)

# Get ranges with high divergence in difference proportion
PropDiff <- ErrorHiFi[ROverlap@from] - ErrorNormal[ROverlap@to]
blnNotNA <- !is.na(PropDiff)
PropDiff[blnNotNA]
plot(PropDiff[!is.na(PropDiff)])

GRHiFi[which.max(ErrorHiFi)]

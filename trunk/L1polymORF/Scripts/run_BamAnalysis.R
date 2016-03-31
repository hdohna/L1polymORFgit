# The script below compares genomic ranges of non-reference L1 peaks from capture and PacBio data

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(ShortRead)
library(csaw)

# Source start script
#source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Summarize results from capture data
CaptureSummary <- BamAnalysis(BamFolder = "home/hzudohna/L1polymORF/Data/NA12878-L15P_NonRef", 
            BamSuffix = "withRG.bam",
            L1Ranges = "home/hzudohna/L1polymORF/Data/NA12878-L15P_L1Ranges.RData",
            L1Length = 6064,
            blnGetSplitReads = T,
            BamFileOriginal = '/home/hzudohna/NA12878-L15P_S1_L001_001.sorted.dedup.mapnonzero.bam')
  
save.image("home/hzudohna/L1polymORF/Data/BamAnalysisOutput.RData")
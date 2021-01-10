# The following script creates a bed file for the non-polymorphic L1 ranges for the reference
# genome hg19 

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(rtracklayer)
library(GenomicRanges)

# Read in table with L1HS on Hg38
L1Table <- read.csv("D:/OneDrive - American University of Beirut/L1polymORF/Data/repeatsHg38_L1HS.csv",
                    as.is = T)
L1Table$ChrNr <- substr(L1Table$genoName, 4, nchar(L1Table$genoName))
L1Table$genoStartMinus1000 <- L1Table$genoStart - 1000
L1Table$genoEndPlus1000    <- L1Table$genoEnd + 1000
L1GR <- makeGRangesFromDataFrame(L1Table, seqnames.field = "genoName",
                                 start.field = "genoStart",
                                 end.field = "genoEnd")
L1GR_left <- makeGRangesFromDataFrame(L1Table, 
                                      seqnames.field = "genoName",
                                      start.field = "genoStartMinus1000",
                                      end.field = "genoStart")

L1GR_right <- makeGRangesFromDataFrame(L1Table, 
                                       seqnames.field = "genoName",
                                       start.field = "genoEnd",
                                       end.field = "genoEndPlus1000")

export.bed(L1GR, "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HSRefRanges_hg38.bed")
export.bed(L1GR_left, "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HSRefRanges_hg38_left1000.bed")
export.bed(L1GR_right, "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HSRefRanges_hg38_right1000.bed")

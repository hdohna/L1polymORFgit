# The following script reads in 1000 genome data on deletions and determines
# which intersect with L1 insertions in the reference genome

# Load packages
library(GenomicRanges)
library(rtracklayer)

# Read in info from 1000 genome deletions and turn them into genomic ranges
DelInfo <- read.delim("/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/Del_chr1_InfoCols",
                      header = F, sep = " ", skip = 1)

DelInfo_GR <- GRanges(seqnames = paste("chr", DelInfo$V1, sep = ""),
        ranges = IRanges(start = DelInfo$V2, end = DelInfo$V2))

# Read in table for reference genome
L1Ref    <- read.csv("/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/repeatsHg19_L1HS.csv")
L1Ref_GR <- makeGRangesFromDataFrame(L1Ref, 
                                     seqnames.field = "genoName",
                                     start.field="genoStart",
                                     end.field="genoEnd")

# Get overlaps between deltions and L1
DelL1RefOverlap <- findOverlaps(DelInfo_GR, L1Ref_GR)

# Calculate deletion size
DelInfo$DelSize <- nchar(as.character(DelInfo$V4)) - nchar(as.character(DelInfo$V5))

# Compare deletion size with L1 insertion size
cbind(DelInfo$DelSize[DelL1RefOverlap@from], width(L1Ref_GR)[DelL1RefOverlap@to])
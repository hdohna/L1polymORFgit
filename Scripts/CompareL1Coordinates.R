# The script below reads files with genomic coordinates of different data sources and
# compares them

##############
# Source prerequisites
##############

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(Rsamtools)
library(rtracklayer)
library(GenomicRanges)

#
################################
#                              #
#  Load and preprocess data    #
#                              #
################################


######
# Read and process repeatmasker table
######

# Read repeat table and subset to get only L1HS rows with fragment size below 
# MaxFragLength
#RepeatTable <- read.delim("D:/L1polymORF/Data/repeatsHg19_L1HS")
# repColNames <- c(paste("X", 1:5, sep = ""), 
#                  c("chromosome", "start", "end", "X9", "strand", "repType", 
#                    "repClass", "repName"), paste("X", 14:17, sep = ""))
# RepeatTable <- read.delim("D:/L1polymORF/Data/rmsk_hg19.txt", header = F,
#                           col.names = repColNames)
# RepeatTable <- RepeatTable[RepeatTable$repName == "L1",]
# nrow(RepeatTable)
# # Create genomic ranges for L1 fragments, match them to distances to get distance
# # to consensus per fragment
# L1RefGR1 <- makeGRangesFromDataFrame(RepeatTable)

# Import bed file with 
L1RefGR2 <- import.bed("D:/L1polymORF/Data/hg19.fa_LINE1.bed")

# Read in table with L1 ranges
L1Table <- read.csv("D:/L1polymORF/Data/L1_repeat_table_Hg19.csv")
L1TableGR <- makeGRangesFromDataFrame(L1Table, seqnames.field = "genoName",
                                      start.field = "genoStart",
                                      end.field = "genoEnd")

# Load L1 references data
load("D:/L1polymORF/Data/L1RefRanges_hg19.Rdata")
load("D:/L1polymORF/Data/BZ_L1Ranges.Rdata")
SuspectL1Ranges

# Create ranges for previously called non-reference insertions
PreviousNonRefL1GR <- GRanges(seqnames = c("chr12", "chr6", "chr7", "chr7", "chrX"),
                              IRanges(start = c(28918987, 157968729, 8019012, 156702046, 71811423),
                                      end = c(28918988, 157968730, 8019013, 156702047, 71811424)))

###############################################
#                                             #
#  Compare different LINE-1 data source       #
#                                             #
###############################################

Dist2Closest(PreviousNonRefL1GR, L1TableGR)
# Get GRanges that are unique to each dataset
blnL1RefGR1inGR2 <- overlapsAny(L1RefGR1, L1RefGR2)
blnL1RefGR2inGR1 <- overlapsAny(L1RefGR2, L1RefGR1)

# Explore unique genomic ranges
UniqueGR1 <- L1RefGR1[!blnL1RefGR1inGR2]
UniqueGR2 <- L1RefGR2[!blnL1RefGR2inGR1]
max(width(UniqueGR1))
max(width(UniqueGR2))

length(UniqueGR2)-length(UniqueGR1)

# Get GRanges that overlap with  
TestRange <- GRanges(seqnames = "chr12", IRanges(start = 28914052, end = 28918987))
idxOverlap2 <- findOverlaps(L1RefGR2, TestRange)
idxOverlap1 <- findOverlaps(L1RefGR1, TestRange)
idxOverlap3 <- findOverlaps(L1TableGR, TestRange)
idxOverlap4 <- findOverlaps(SuspectL1Ranges, TestRange)
findOverlaps(L1GRanges, TestRange)
length(L1RefGR2) - length(L1RefGR1)
L1RefGR1[idxOverlap1@from]
RepeatTable[idxOverlap1@from,]
sum(!blnL1RefGR1inGR2)
sum(!blnL1RefGR2inGR1)
table(L1RefGR2@elementMetadata@listData$name[!blnL1RefGR2inGR1])
hist(width(L1RefGR2[!blnL1RefGR2inGR1]))
sum(!blnL1RefGR2inGR1)

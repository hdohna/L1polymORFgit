##############################################
#
# General description:
#
#   The following script reads in genomic hg19 ranges of PacBio reads that map
#   onto L1HS consensus sequence and determines which of these ranges overlap
#   with reference L1.

# Input:
#
#     BamFileToBeFiltered: path to bam file that has reads mapped to hg19
#     BamFilter: path to bam file that contains reads mapped to L1HS

# Output:
#   
#    : ...

##############################################

######                                      
# Source packages and set parameters  
######

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')
#source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(ShortRead)
library(Rsamtools)
library(csaw)
library(GenomicRanges)

# Files and folders
L1TableFileName <- "/home/hzudohna/L1polymORF/Data/L1_repeat_table_Hg19.csv"
OutFile         <- "/home/hzudohna/L1polymORF/Data/AnalyzedPacBioL1Ranges.RData"

#L1TableFileName <- "D:/L1polymORF/Data/L1_repeat_table_Hg19.csv"

# Peak calling parameters
MinMaxCover <- 5    # minimum maximum coverage to be called a peak 
MinGap      <- 10
MinDist2L1  <- 3*10^4 # minimum distance to L1 to be called a peak 
PacBioWindow <- 1000

# Load genomic ranges of reads mapped to L1HS
load(file = "/home/hzudohna/L1polymORF/Data/PacBioL1Ranges.RData")
#load(file = "D:/L1polymORF/Data/PacBioL1Ranges.RData")

#######################################
#                                     #
#    Read in data and                 #
#   determine coverage                #
#                                     #
#######################################

# Read in table with L1 ranges
L1Table <- read.csv(L1TableFileName)

# Create GRanges objects with L1 Seqences
L1IRanges <- IRanges(start = L1Table$genoStart,
                     end = L1Table$genoEnd)
L1GRanges <- GRanges(seqnames = L1Table$genoName, ranges = L1IRanges,
                     strand = L1Table$strand)

cat("*******   Turning BAM files into GRanges ...   *******\n")

# Read coverage per chromosome
CoverList <- lapply(1:length(ReadsPerChromPacBioL1), function(i){
  
  cat("Calculating coverage for chromosome", i, "\n")
  Reads   <- ReadsPerChromPacBioL1[[i]]
  ReadCov <- coverage(Reads)
})

#######################################
#                                     #
#    Determine 'islands' and          #
#        overlap with L1              #
#                                     #
#######################################

# Determine separate islands with continuous read coverage and turn islands 
# into genomic ranges
IslandList <- lapply(CoverList, function(x){
  Islands <- slice(x, lower = 1)
})
length(IslandList[[1]]@listData$chr1@subject@values)
Chroms <- paste("chr", c(1:22, "X", "Y"), sep = "")
IslandGRanges <- lapply(1:length(IslandList), function(i){
  GRanges(seqnames = Chroms[i], 
          ranges = IslandList[[i]]@listData[[1]]@ranges,
          coverTotal = viewSums(IslandList[[i]])[[1]],
          coverMax   = viewMaxs(IslandList[[i]])[[1]],
          coverMaxPos   = viewWhichMaxs(IslandList[[i]])[[1]])
})
IslandGRanges <- GRangesList(IslandGRanges)
IslandGRanges <- unlist(IslandGRanges)
cat(length(IslandGRanges), "distinct peaks\n\n")

# Merge ranges that are less than MinGap bp apart
IslGRanges_reduced <- reduce(IslandGRanges, min.gapwidth = MinGap,
                             with.revmap = T)
cat(length(IslGRanges_reduced), "distinct peaks after merging peaks that\n")
cat("are less than", MinGap, "apart\n\n")

# Find overlaps between islands and L1HS ranges
blnOverlapIslands_All <- overlapsAny(IslGRanges_reduced, L1GRanges)

# Get ranges of suspected L1s
maxCoverOriginal    <- IslandGRanges@elementMetadata@listData$coverMax
maxCoverPosOriginal <- IslandGRanges@elementMetadata@listData$coverMaxPos
maxCover <- sapply(IslGRanges_reduced@elementMetadata@listData$revmap, 
                   function(x) max(maxCoverOriginal[x]))
maxCoverPos <- sapply(IslGRanges_reduced@elementMetadata@listData$revmap, 
                      function(x) maxCoverPosOriginal[x[which.max(maxCoverOriginal[x])]])
idxSuspectL1Ranges <- which(maxCover > MinMaxCover & (!blnOverlapIslands_All))
SuspectL1Ranges    <- IslGRanges_reduced[idxSuspectL1Ranges]
cat("\n", length(idxSuspectL1Ranges), "peaks have maximum coverage of at least",
    MinMaxCover, "and do not overlap with reference L1\n")

# Remove ranges of suspected L1s that are too close
DistToNearestL1    <- nearest(SuspectL1Ranges, L1GRanges)
idxSuspectL1Ranges <- idxSuspectL1Ranges[DistToNearestL1 >= MinDist2L1]
SuspectL1Ranges    <- IslGRanges_reduced[idxSuspectL1Ranges]
maxCoverPos_SuspL1Ranges <- maxCoverPos[idxSuspectL1Ranges]
cat(length(idxSuspectL1Ranges), "of the above peaks are at least", MinDist2L1,
    "bp from nearest reference L1\n\n")

#######################################################
#                                                     #
#    Save results                                     #
#                                                     #
#######################################################

# Save results
cat("*******  Saving results ...   *******\n")
save.image(file = OutFile)


##############################################
#
# General description:
#
#   The following script reads a bam file, finds peaks that do not overlap
#   with reference L1s 

# Input:
#
#     BamFile: path to file that contains mapped reads
#     L1TableFileName: path to file that contains L1HS ranges in a table

# Output:
#   
#    : ...

##############################################

########################################
#                                      #
#  Source packages and set parameters  #
#                                      #
########################################

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')
#source('D:/L1polymORF/Scripts/_Start_L1polymORF.R')

# Load packages
library(ShortRead)
library(csaw)
library(chipseq)
library(BSgenome.Hsapiens.UCSC.hg38)

# Get all ranges of reads for per chromosome
Chroms       <- paste('chr', c(1:22, "X", "Y"), sep = "")

# Files and folders
BamFile         <- "/home/hzudohna/BoData/NA12878-L1HS_S1_L001.dedup.unique.sorted.bam"
L1TableFileName <- "/home/hzudohna/L1polymORF/Data/L1_repeat_table.csv"
ChainFilePath   <- "/home/hzudohna/L1polymORF/Data/hg38ToHg19.over.chain"
OutResults      <- '/home/hzudohna/L1polymORF/Data/NonRefL1Ranges.Rdata'

# BamFile         <- "D:/L1polymORF/Data/NA12878-L1HS_S1_L001.dedup.unique.sorted.bam"
# L1TableFileName <- "D:/L1polymORF/Data/L1_repeat_table.csv"
# ChainFilePath   <- "D:/L1polymORF/Data/hg38ToHg19.over.chain"
# OutResults      <- 'D:/L1polymORF/Data/NonRefL1Ranges.Rdata'


# Peak calling parameters
MinMaxCover <- 5    # minimum maximum coverage to be called a peak 
MinGap      <- 6000
MinDist2L1  <- 3*10^4 # minimum distance to L1 to be called a peak 

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
CoverList <- lapply(Chroms, function(Chrom){
  cat("Reading reads for chromosome", Chrom, "\n")
  ChromLength <- length(BSgenome.Hsapiens.UCSC.hg38[[Chrom]])
  R1 <- GRanges(seqnames = Chrom, ranges = IRanges(start = 1, 
                                                   end = ChromLength))
  Reads   <- extractReads(bam.file = BamFile, region = R1)
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
IslandGRanges <- lapply(1:length(IslandList), function(i){
  GRanges(seqnames = Chroms[i], 
          ranges = IslandList[[i]]@listData[[1]]@ranges,
          coverTotal = viewSums(IslandList[[i]])[[1]],
          coverMax   = viewMaxs(IslandList[[i]])[[1]],
          coverMaxPos   = viewWhichMaxs(IslandList[[i]])[[1]])
})
IslandGRanges <- GRangesList(IslandGRanges)
IslandGRanges <- unlist(IslandGRanges)

# Merge ranges that are less than MinGap bp apart
IslGRanges_reduced <- reduce(IslandGRanges, min.gapwidth = MinGap,
                        with.revmap = T)

# Find overlaps between islands and L1HS ranges
blnOverlapIslands_All <- overlapsAny(IslGRanges_reduced, L1GRanges)

# Get ranges of suspected L1s
maxCoverOriginal    <- IslandGRanges@elementMetadata@listData$coverMax
maxCoverPosOriginal <- IslandGRanges@elementMetadata@listData$coverMaxPos
maxCover <- sapply(IslGRanges_reduced@elementMetadata@listData$revmap, 
                   function(x) max(maxCoverOriginal[x]))
maxCoverPos <- sapply(IslGRanges_reduced@elementMetadata@listData$revmap, 
   function(x) maxCoverPosOriginal[which.max(maxCoverOriginal[x])])
idxSuspectL1Ranges <- which(maxCover > MinMaxCover & (!blnOverlapIslands_All))
SuspectL1Ranges    <- IslGRanges_reduced[idxSuspectL1Ranges]
cat("\n", length(idxSuspectL1Ranges), "peaks have maximum coverage of at least",
    MinMaxCover, "and do not overlap with reference L1\n")

# Remove ranges of suspected L1s that are too close
DistToNearestL1    <- nearest(SuspectL1Ranges, L1GRanges)
idxSuspectL1Ranges <- idxSuspectL1Ranges[DistToNearestL1 >= MinDist2L1]
SuspectL1Ranges    <- IslGRanges_reduced[idxSuspectL1Ranges]
maxCoverPos_SuspL1Ranges <- maxCoverPos[idxSuspectL1Ranges]
cat(length(idxSuspectL1Ranges), "of the above peaksare at least", MinDist2L1,
    "bp from nearest reference L1\n\n")

# Match coordinates to hg19
SuspL1PeakRanges <- GRanges(seqnames = seqnames(SuspectL1Ranges),
                            ranges = IRanges(start = maxCoverPos_SuspL1Ranges,
                                             end = maxCoverPos_SuspL1Ranges))
SuspectL1Ranges19 <- liftOver(SuspL1PeakRanges, 
                             chain = import.chain(ChainFilePath))
NrMapped <- sapply(SuspectL1Ranges19, length)

# Retain only the coordinates that are uniquely mapped
SuspectL1Ranges19Mapped <- unlist(SuspectL1Ranges19[NrMapped == 1])
idxMapped2hg19 <- idxSuspectL1Ranges[NrMapped == 1]

#######################################################
#                                                     #
#    Save results                                     #
#                                                     #
#######################################################

# Save results
cat("*******  Saving results ...   *******\n")
save.image(file = OutResults)

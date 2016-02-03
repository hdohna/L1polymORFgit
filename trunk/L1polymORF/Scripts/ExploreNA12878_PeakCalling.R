##############################################
#
# General description:
#
#   The following script reads a bam for 12878 and explores peak
#   overlap 

# Input:
#
#     BamFile: path to file that contains mapped reads
#     L1HSTableFileName: path to file that contains L1HS ranges in a table

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
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(pROC)
library(ShortRead)
library(csaw)
library(chipseq)
library(BSgenome.Hsapiens.UCSC.hg19)

# Get all ranges of reads for per chromosome
Chromosomes <- paste("chr", 1:22, sep = "")
MinCoverage <- 5
IslRange    <- 10000
CalcNew_IslandPerCoverage <- T

# Files and folders
BamFile            <- "D:/L1polymORF/Data/NA12878-L1HS_S1_L001.dedup.unique.sorted.bam"
L1HSTableFileName  <- "D:/L1polymORF/Data/L1HS_repeat_table.csv"

#######################################
#                                     #
#    Read in data                     #
#                                     #
#######################################

# Read in table with L1 ranges
L1HSTable <- read.csv(L1HSTableFileName)

# Make some corrections to create a proper GRanges object with L1 Seqences
L1IRanges <- IRanges(start = L1HSTable$genoStart,
                     end = L1HSTable$genoEnd)
L1GRanges <- GRanges(seqnames = L1HSTable$genoName, ranges = L1IRanges,
                     strand = L1HSTable$strand)
L1GRanges <- L1GRanges[width(L1GRanges) >= 6000]

cat("*******   Turning BAM files into GRanges ...   *******\n")

# Read coverage per chromosome
Chroms       <- paste('chr', c(1:22, "X", "Y"), sep = "")
CoverList <- lapply(Chroms, function(Chrom){
  cat("Reading reads for chromosome", Chrom, "\n")
  ChromLength <- length(BSgenome.Hsapiens.UCSC.hg19[[Chrom]])
  R1 <- GRanges(seqnames = Chrom, ranges = IRanges(start = 1, end = ChromLength))
  Reads <- extractReads(bam.file = BamFile , region = R1)
  ReadCov <- coverage(Reads)
})

#######################################
#                                     #
#    Determine 'islands' and          #
#        overlap with L1              #
#                                     #
#######################################

# Determine separate islands with continuous read coverage and turn islands into
# genomic ranges
IslandList <- lapply(CoverList, function(x){
  Islands <- slice(x, lower = 1)
})
IslandGRanges <- lapply(1:length(IslandList), function(i){
  GRanges(seqnames = Chroms[i], 
          ranges = IslandList[[i]]@listData[[1]]@ranges,
          coverTotal = viewSums(IslandList[[i]])[[1]],
          coverMax   = viewMaxs(IslandList[[i]])[[1]])
})
IslandGRanges <- GRangesList(IslandGRanges)
IslandGRanges <- unlist(IslandGRanges)
RangeHist1    <- hist(width(IslandGRanges), breaks = seq(0, 25000, 100), 
                   plot = F)
IslandGRanges@elementMetadata@listData$coverTotal

# Get the coverage sum and maximum
SumList <- lapply(IslandList, function(x) viewSums(x))
MaxList <- lapply(IslandList, function(x) viewMaxs(x))
names(SumList[[1]])

# Find overlaps between islands and L1HS ranges
IslandL1Overlaps  <- findOverlaps(IslandGRanges, L1GRanges)
blnOverlapIslands <- 1:length(IslandGRanges) %in% IslandL1Overlaps@queryHits
blnOverlapL1      <- 1:length(L1GRanges) %in% IslandL1Overlaps@subjectHits

# Subset to get islands overlapping with L1HS
IslandGRangesSubset <- IslandGRanges[blnOverlapIslands]
RangeHist2 <- hist(width(IslandGRangesSubset), breaks = seq(0, 25000, 100), 
                   plot = F)
plot(RangeHist1$mids, RangeHist1$density, type = "l", ylab = "Density", 
     xlab = "Width [bp]", ylim = c(0, 0.0005))
lines(RangeHist2$mids, RangeHist2$density, col = "red")

# Get ranges around each L1Hs and load intersecting reads
LargeL1Ranges <- resize(L1GRanges,10^4, fix = "center")
ReadsPerL1 <- lapply(LargeL1Ranges, function(x) extractReads(bam.file = BamFile , x))

#######################################
#                                     #
#    Explore peak calling             #
#                                     #
#######################################

# Determine receiver-operating curves based on total and maximum coverage 
cat("*******   Calculating ROC curves ...   *******\n")
# rocTotal <- pROC::roc(response = blnOverlapIslands,
#                 predictor = IslandGRanges@elementMetadata@listData$coverTotal)
# rocMax   <- pROC::roc(response = blnOverlapIslands,
#                 predictor = IslandGRanges@elementMetadata@listData$coverMax)
# 
# # Plot ROCs
# Cols <- rainbow(2)
# plot(1 - rocTotal$specificities, rocTotal$sensitivities,
#      xlab = "1 - specificity", ylab = "Sensitivity", col = Cols[1])
# points(1 - rocMax$specificities, rocMax$sensitivities, col = Cols[2])
# legend("bottomright", legend = c("Total coverage", "Maximum coverage"),
#        col = Cols, lty = c(1,1))
# CreateDisplayPdf('D:/L1polymORF/Figures/PeakCallingROC.pdf')
# 
# rocMax$thresholds[rocMax$specificities > 0.9 & rocMax$sensitivities > 0.6]
hist(IslandGRanges@elementMetadata@listData$coverMax,
     breaks = 0:2500, xlim = c(0, 300), ylim  = c(0, 1000))

IslandGRanges@elementMetadata@listData$coverMax > 100

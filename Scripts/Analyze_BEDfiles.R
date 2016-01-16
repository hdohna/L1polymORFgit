##############################################
#
# General description:
#
#   The following script reads bam files of reads coming from chip-seq with
#   capture oligos containing L1HS sequences

# Input:
#
#    hg38.fa_L1HS_6kb.bed: bed file with full-length L1
#    NA12878_L1_capt_IDTx_peaks.bed: bed file with peaks on Chip-seq on NA12878

# Output:
#   
#    : ...

##############################################

# Source start script
#source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(ShortRead)
library(csaw)
library(chipseq)
library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(VennDiagram)


# Get all ranges of reads for per chromosome
Chromosomes <- paste("chr", 1:22, sep = "")
MinCoverage <- 5
IslRange    <- 1000
CalcNew_IslandPerCoverage <- T

#######################################
#                                     #
#    Turn BAM files into GRanges      #
#                                     #
#######################################

cat("*******   Turning BED files into GRanges ...   *******\n")
L1Ref        <- import.bed(con = "D:/L1polymORF/Data/hg38.fa_L1HS_6kb.bed") 
NA12878Peaks <- import.bed(con = "D:/L1polymORF/Data/NA12878_L1_capt_IDTx_peaks.bed") 

# Reduce peaks by joining ranges that are separated by less than 7000 bp
NA12878PeaksReduced <- reduce(NA12878Peaks, min.gapwidth = 5000)

# Relationship between window around a peak and number of peaks intersecting with
WindowSizes <- seq(1000, 10000, 1000)
NrIntersect <- sapply(WindowSizes, function(x) {
  NewPeaks <- NA12878Peaks
  ranges(NewPeaks) <- resize(ranges(NA12878Peaks), width = x, fix = "center") 
  sum(L1Ref %over% NewPeaks )
})
plot(WindowSizes, NrIntersect)

# Relationship between window around a peak and number of peaks intersecting with
hist(width(NA12878Peaks), xlab = "Peak width [bp]", main = "")
dev.copy2pdf(file = "D:/L1polymORF/Figures/PeakWidth.pdf")
max(width(NA12878PeaksReduced))
BinWidth <- 1500
PeakWidths <- seq(0, 9000, BinWidth)
StandardPeaks <- NA12878PeaksReduced
ranges(StandardPeaks) <- resize(ranges(StandardPeaks), 2000, fix = "center") 
PropIntersect <- sapply(2:length(PeakWidths), function(i) {
  WidthMask <- width(NA12878PeaksReduced) >= PeakWidths[i - 1] &
    width(NA12878PeaksReduced) < PeakWidths[i]
  NewPeaks <- StandardPeaks[WidthMask]
  sum(NewPeaks %over% L1Ref ) / sum(WidthMask)
})

plot(PeakWidths[-1] - 0.5 * BinWidth, PropIntersect, 
     xlab = "Proportion of peaks intersecting",
     ylab = "Peak width")
dev.copy2pdf(file = "D:/L1polymORF/Figures/PeakIntersectionVsLength.pdf")


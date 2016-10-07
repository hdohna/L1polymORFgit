# The following script test for each L1 in the catalog whether it is present
# in a genome using PacBio data from a capture experiment

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(csaw)
library(ShortRead)
library(Rsamtools)

# Specify flank width
Fwidth <- 50

# Specify path to PacBio bam file (reads alinged to 
# L1CatalogueWithFlank_Sat_May_07_15-15-31_2016.fa)
BamFilePath      <- "D:/L1polymORF/Data/NA12878_L1PCR_hg19withL1.sorted.bam"
list.files("D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_ROI_catl110kb/")

# Read file with PCR results
L1_PCR_results <- read.delim("D:/L1polymORF/Data/L1_PCR_results.txt", header = F,
                             as.is = T, skip = 2, sep = " ", 
                             col.names = c("Name", "Well", "Product"))

# Read in table with known L1 
L1Catalog   <- read.csv("D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv",
                      as.is = T)
blnL1Mapped  <- !is.na(L1Catalog$start_HG38)
blnL1Allele1 <- L1Catalog$Allele == 1
L1CatalogMapped <- L1Catalog[blnL1Mapped & blnL1Allele1,]

# Read in ranges 
L1RangesDF   <- read.table("D:/L1polymORF/Data/L1Ranges_hg19_withNonRefL1", as.is = T,
                           header = T)

# Turn L1 ranges into a GRanges object
L1Ranges <- makeGRangesFromDataFrame(L1RangesDF)
L1RangesLeft  <- flank(L1Ranges, width = 200,  ignore.strand = T)
L1RangesRight <- flank(L1Ranges, width = 200, start = F,
                       ignore.strand = T)
L1RangesLeft  <- shift(L1RangesLeft, -50)
L1RangesRight <- shift(L1RangesRight, 50)

# Loop through flanking ranges elements and extract reads
cat("Extracting reads from bam file ... \n")
NrOverlap_Left  <- rep(NA, length(L1Ranges))
NrOverlap_Right <- rep(NA, length(L1Ranges))
for(i in 1:length(L1Ranges)){
  cat("Extracting reads for L1",i,"\n")
  RLeft  <- extractReads(BamFilePath, L1RangesLeft[i])
  RRight <- extractReads(BamFilePath, L1RangesRight[i])
  NrOverlap_Left[i]  <-  length(RLeft)
  NrOverlap_Right[i] <- length(RRight)
}
# MaxRangeList <- lapply(seq_along(L1withFlankRanges), function(i){
#   L1Range <- L1withFlankRanges[i]
#   L1Name <- as.vector(seqnames(L1withFlankRanges))[i]
#   Reads <- extractReads(bam.file = BamFilePath, region = L1Range)
#   ReadCov <- coverage(Reads)
#   Islands <- slice(ReadCov, lower = 1)
#   GRs <- GRanges(seqnames = L1Name, 
#                  ranges = Islands@listData[[1]]@ranges,
#                  coverTotal = viewSums(Islands)[[1]],
#                  coverMax   = viewMaxs(Islands)[[1]],
#                  coverMaxPos   = viewWhichMaxs(Islands)[[1]])
#   idxMaxTot <- which.max(GRs@elementMetadata@listData$coverTotal)
#   GRs[idxMaxTot]
# })
# blnGRanges <- sapply(ReadRanges, function(x)is(x)[1] == "GRanges")
# ReadRanges <- GRangesList(ReadRanges[blnGRanges])
# ReadRanges <- unlist(ReadRanges)
# 
# # Define flank ranges for getting reads intersecting with
# LeftFlankRanges   <- GRanges(seqnames = colnames(L1StartEnd),
#                              ranges = IRanges(start = L1StartEnd["Start", ] - 300, 
#                                               end = L1StartEnd["Start", ] - 5))
# RightFlankRanges  <- GRanges(seqnames = colnames(L1StartEnd),
#                              ranges = IRanges(start = L1StartEnd["End", ] + 5, 
#                                               end = L1StartEnd["End", ] + 300))
# 
# # Count overlaps for left and right 
# NrOverlap_Left  <- countOverlaps(LeftFlankRanges, ReadRanges, minoverlap = 1)
# NrOverlap_Right <- countOverlaps(RightFlankRanges, ReadRanges, minoverlap = 1)
# ReadRanges[NrOverlap_Left]
# 

# Get accession number from file with primer results and match to catalog
L1_PCR_results$AccNr <- sapply(L1_PCR_results$Name, 
                               function(x) strsplit(x, "_")[[1]][1])
AccMatch <- match(L1_PCR_results$AccNr, L1RangesDF$name)
NrOverlap_Left_matched  <- NrOverlap_Left[AccMatch]
NrOverlap_Right_matched <- NrOverlap_Right[AccMatch]

# Get all L1 ranges that were not part of PCR
blnNoPCR <- !L1RangesDF$name %in% L1_PCR_results$AccNr
NrOverlap_NoPCR <- c(NrOverlap_Left[blnNoPCR], NrOverlap_Right[blnNoPCR])
DensNoPCR <- density(NrOverlap_NoPCR, from = 0)

# Get for each L1_PCR_results the number of matching reads
NcharAcc <- nchar(L1_PCR_results$AccNr)
L1_PCR_results$Side <- substr(L1_PCR_results$Name, NcharAcc + 3, NcharAcc + 3)
NrOverlap_matched <- sapply(1:nrow(L1_PCR_results), function(x){
  if (L1_PCR_results$Side[x] == "1"){
    NrOverlap_Left_matched[x]
  } else {
    NrOverlap_Right_matched[x]
  }
})
L1_PCR_results$NrReadsInFlank <- NrOverlap_matched

# Analyze relationship between number of overlapping reads and PCR product
table(L1_PCR_results$Product, NrOverlap_matched > 0)
table(L1_PCR_results$Product, NrOverlap_matched > 100)
boxplot(NrOverlap_matched ~ L1_PCR_results$Product)

# Determine ratios of 1- cumulative densities
idxNotNA <- which(!is.na(NrOverlap_matched))
idxSimilar <- sapply(NrOverlap_matched[idxNotNA], 
                     function(x) which.min(abs(x - NrOverlap_NoPCR)))
F_PCR <- sapply(NrOverlap_matched[idxNotNA], function(x) {
  sum(x <= NrOverlap_matched[idxNotNA])}) / length(idxNotNA)
F_NoPCR <- sapply(NrOverlap_NoPCR, function(x) {
  sum(x <= NrOverlap_NoPCR)}) / length(NrOverlap_NoPCR)
ProbReal <- 1 - F_NoPCR[idxSimilar] / F_PCR 
mean(NrOverlap_NoPCR > 0)
ProbReal[NrOverlap_matched[idxNotNA] > max(NrOverlap_NoPCR)] <- 1

# Add posterior probability and save PCR result file
L1_PCR_results$ProbPresent <- NA
L1_PCR_results$ProbPresent[idxNotNA] <- ProbReal
write.csv(L1_PCR_results, "D:/L1polymORF/Data/NA12878_L1PCR_results.csv")
plot(L1_PCR_results$NrReadsInFlank, L1_PCR_results$ProbPresent)


# # Get highest peak per L1
# i <- 1
# MaxRangeList <- lapply(seq_along(L1withFlankRanges), function(i){
#   L1Range <- L1withFlankRanges[i]
#   L1Name <- as.vector(seqnames(L1withFlankRanges))[i]
#   Reads <- extractReads(bam.file = BamFilePath, region = L1Range)
#   ReadCov <- coverage(Reads)
#   Islands <- slice(ReadCov, lower = 1)
#   GRs <- GRanges(seqnames = L1Name, 
#                  ranges = Islands@listData[[1]]@ranges,
#                  coverTotal = viewSums(Islands)[[1]],
#                  coverMax   = viewMaxs(Islands)[[1]],
#                  coverMaxPos   = viewWhichMaxs(Islands)[[1]])
#   idxMaxTot <- which.max(GRs@elementMetadata@listData$coverTotal)
#   GRs[idxMaxTot]
# })
# MaxRangeList <- GRangesList(MaxRangeList)
# MaxRanges    <- unlist(MaxRangeList)
# 
# # Determine whether highest peak intersect with left, right or no end
# all(as.vector(seqnames(L1Ranges)) == as.vector(seqnames(L1withFlankRanges)))
# PeakOverlap <- sapply(seq_along(MaxRanges), function(i) {
#   RangeEnds <- c(start(L1Ranges[i]), end(L1Ranges[i]))
#   L1Name <- as.vector(seqnames(L1withFlankRanges))[i]
#   CompareRanges <- GRanges(seqnames = L1Name, 
#                            ranges = IRanges(start = RangeEnds, end = RangeEnds))
#   idxOverlap <- which(overlapsAny(CompareRanges, MaxRanges[i]))
#   if (length(idxOverlap) == 0) idxOverlap <- 3
#   idxOverlap
# })
# 
# PeakCorrect <- L1_PCR_results$Side == as.character(PeakOverlap[AccMatch])
# table(L1_PCR_results$Product, PeakCorrect)


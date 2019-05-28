# Load packages
library(GenomicRanges)

load("D:/L1polymORF/Data/L1_TajimaD.RData")

# Create genomic ranges 
TajimaD_GR <- GRanges()
TajimaD_reordered <- data.frame()

# Loop over chromosomes, order 
for (chr in unique(TajimaData$CHROM)){
  TDataSubset <- TajimaData[TajimaData$CHROM == chr,]
  startOrder  <- order(TDataSubset$BIN_START)
  TDataSubset <- TDataSubset[startOrder, ]
  New_GR <- GRanges(seqnames = chr, 
                        ranges = IRanges(start = TDataSubset$BIN_START[-nrow(TDataSubset)],
                                         end  = TDataSubset$BIN_START[-1]))  
  TDataSubset    <- TDataSubset[-1, ]
  
  # Retain only GRs and TData rows where New_GR have proper width
  blnProperWidth <- width(New_GR) == StepSize + 1
  New_GR         <- New_GR[blnProperWidth]
  TDataSubset    <- TDataSubset[blnProperWidth, ]
  
  # Append new data
  TajimaD_reordered <- rbind(TajimaD_reordered, TDataSubset)
  TajimaD_GR <- c(TajimaD_GR, New_GR)
}
TajimaD_GR <- as(TajimaD_GR, "GRanges")

# Find ranges with only one L1
TajDL1_OLcount <- countOverlaps(TajimaD_GR, L1GRanges)
TajimaD_GR1 <- TajimaD_GR[TajDL1_OLcount == 1]

# Get indices of L1 overlapping with TajimaD_GR1
idxL1_OL <- which(overlapsAny(L1GRanges, TajimaD_GR1))

# Get indices of L1s with non-overlapping neighborhoods
idxUniqueL1Neighbor <- idxL1_OL[countOverlaps(L1Neighborhoods[idxL1_OL], 
                                           L1Neighborhoods[idxL1_OL]) == 1]

# Find overlaps
OL_L1TD <- findOverlaps(L1Neighborhoods[idxUniqueL1Neighbor], TajimaD_GR)
TDList <- lapply(unique(OL_L1TD@from), function(x) {
  idxNeibhborL1 <- OL_L1TD@to[OL_L1TD@from == x]
  TDs <- TajimaD_reordered$TajimaD[idxNeibhborL1]
  TDs - mean(TDs[c(1, length(TDs))])
})
table(sapply(TDList, length))

# Create a matrix of Tajima's D values
NrInt <- FlankSize %/% StepSize
TDMat <- sapply(TDList, function(x){
  if (length(x) < NrInt){
    rep(NA, NrInt)
  } else {
    x[1:NrInt]
  }
})

###########
# Plot mean Tajima's D per Fragments and full-length L1
###########

# Indices of fragment and full-length L1
idxFull <- which(width(L1GRanges[idxUniqueL1Neighbor][unique(OL_L1TD@from)]) >= 6000)
idxFragm <- setdiff(1:length(unique(OL_L1TD@from)), idxFull)

# Mean Tajima's D per fragment and full-length L1
MeanTD_Full  <- rowMeans(TDMat[,idxFull], na.rm = T)
MeanTD_Fragm <- rowMeans(TDMat[,idxFragm], na.rm = T)
MedTD_Full  <- apply(TDMat[,idxFull], 1, FUN = function(x) median(x, na.rm = T))
MedTD_Fragm <- apply(TDMat[,idxFragm], 1, FUN = function(x) median(x, na.rm = T))

# Plot average Tajima's D for fragment and full-length L1 
plot(MeanTD_Full, type = "l", col = "red", xlab = "Genomic position (relative to L1)",
     ylab = "Mean Tajima's D", xaxt = "n")
axis(1, at = seq(0, 50, 10), labels = seq(0, FlankSize, 10*StepSize) - 0.5*FlankSize)
lines(MeanTD_Fragm, col = "blue")

# Plot average Tajima's D for fragment and full-length L1 
plot(MedTD_Full, type = "l", col = "red", xlab = "Genomic position (relative to L1)",
     ylab = "Median Tajima's D", xaxt = "n")
axis(1, at = seq(0, 50, 10), labels = seq(0, FlankSize, 10*StepSize) - 0.5*FlankSize)
lines(MedTD_Fragm, col = "blue")


# Test for different Tajima's D between full-length and fragment at the midpoint
MidPt  <- NrInt %/% 2
t.test(TDMat[MidPt,idxFull], TDMat[MidPt,idxFragm])

# Trend in mean Tajima's D
LMTrend_Fragm <- lm(MeanTD_Fragm ~ c(1:length(MeanTD_Fragm)))
summary(LMTrend_Fragm)
LMTrend_Full <- lm(MeanTD_Full[-MidPt] ~ c(1:length(MeanTD_Full))[-MidPt])
summary(LMTrend_Full)

# Trend in median Tajima's D
LMMedTrend_Fragm <- lm(MedTD_Fragm ~ c(1:length(MedTD_Fragm)))
summary(LMMedTrend_Fragm)
LMMedTrend_Full <- lm(MedTD_Full[-MidPt] ~ c(1:length(MedTD_Full))[-MidPt])
summary(LMMedTrend_Full)

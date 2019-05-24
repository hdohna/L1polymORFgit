# Load packages
library(GenomicRanges)

load("D:/L1polymORF/Data/L1_TajimaD.RData")

# Create genomic ranges 
TajimaD_GR <- GRanges()
TajimaD_reordered <- data.frame()
for (chr in unique(TajimaData$CHROM)){
  TDataSubset <- TajimaData[TajimaData$CHROM == chr,]
  startOrder <- order(TDataSubset$BIN_START)
  TDataSubset <- TDataSubset[startOrder, ]
  New_GR <- GRanges(seqnames = chr, 
                        ranges = IRanges(start = TDataSubset$BIN_START[-nrow(TDataSubset)],
                                         end  = TDataSubset$BIN_START[-1]))  
  TajimaD_GR <- c(TajimaD_GR, New_GR)
  TajimaD_reordered <- rbind(TajimaD_reordered, TDataSubset[-1, ])
}
TajimaD_GR <- as(TajimaD_GR, "GRanges")
blnProperWidth <- width(TajimaD_GR) == StepSize + 1
TajimaD_GR <- TajimaD_GR[blnProperWidth]
TajimaD_reordered <- TajimaD_reordered[blnProperWidth, ]

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
  blnL1 <- OL_L1TD@from == x
  TDs <- TajimaD_reordered$TajimaD[blnL1]
  TDs - mean(TDs[c(1, length(TDs))])
})
table(sapply(TDList, length))
plot(TDList[[3]], type = "l", col = rgb(0, 0, 0, alpha = 0.1),
     ylim = c(-2.3, 2))
for (i in 1:500){
  lines(TDList[[i]], col = rgb(0, 0, 0, alpha = 0.01))
}

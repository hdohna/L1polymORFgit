# The following script reads data of Tajima's D around reference L1 (created 
# in script TajimasDAroundL1.R) and plots mean and median Tajima's D around
# fragment and full-length L1

# Load packages
library(GenomicRanges)

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Read in data (created in script TajimasDAroundL1.R)
load("D:/L1polymORF/Data/L1_TajimaData.RData")

# List of data objects with Tajima's D
TDList_AllPops <- ProcessTajimaDData(L1GRanges = L1GRanges,
                                     StepSize = StepSize,
                                     FlankSize = FlankSize)


# Add a column with NA for intoerval
TajimaData$TajimaD_modified <- TajimaData$TajimaD
TajimaData$TajimaD_modified[TajimaData$N_SNPS == 0] <- NA

# Intialize genomic ranges and matching dataframe with Tajima's D
TajimaD_GR <- GRanges()
TajimaD_reordered <- data.frame()

# Loop over chromosomes, order positions per chromosome
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
TajimaD_GR1    <- TajimaD_GR[TajDL1_OLcount == 1]

# Get indices of L1 overlapping with TajimaD_GR1
idxL1_OL <- which(overlapsAny(L1GRanges, TajimaD_GR1))

# Get indices of L1s with non-overlapping neighborhoods
idxUniqueL1Neighbor <- idxL1_OL[countOverlaps(L1Neighborhoods[idxL1_OL], 
                                           L1Neighborhoods[idxL1_OL]) == 1]

# Find overlaps
OL_L1TD <- findOverlaps(L1Neighborhoods[idxUniqueL1Neighbor], TajimaD_GR)

# Create a list of Tajima's D values around each 
TDList <- lapply(unique(OL_L1TD@from), function(x) {
  idxNeibhborL1 <- OL_L1TD@to[OL_L1TD@from == x]
  TDs <- TajimaD_reordered$TajimaD_modified[idxNeibhborL1]
  EdgeMean <- mean(TDs[c(1, length(TDs))], na.rm = T)
  if (is.na(EdgeMean)) {EdgeMean <- 0}
  TDs - EdgeMean
})

# Create a matrix of Tajima's D values
NrInt <- FlankSize %/% StepSize
MidPt <- NrInt %/% 2
DiffFromCenter <- 1:NrInt - MidPt
AbsDiff <- abs(DiffFromCenter)
sapply(TDList, length)
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
MedTD_Full   <- apply(TDMat[,idxFull], 1, FUN = function(x) median(x, na.rm = T))
MedTD_Fragm  <- apply(TDMat[,idxFragm], 1, FUN = function(x) median(x, na.rm = T))

# Plot average Tajima's D for fragment and full-length L1 
plot(MeanTD_Full, type = "l", col = "red", xlab = "Genomic position (relative to L1)",
     ylab = "Mean Tajima's D", xaxt = "n")
axis(1, at = seq(0, 50, 10), labels = seq(0, FlankSize, 10*StepSize) - 0.5*FlankSize)
lines(MeanTD_Fragm, col = "blue")

# Plot median Tajima's D for fragment and full-length L1 
plot(MedTD_Full, col = "red", xlab = "Genomic position (relative to L1)",
     ylab = "Median Tajima's D", xaxt = "n")
axis(1, at = seq(0, 50, 10), labels = seq(0, FlankSize, 10*StepSize) - 0.5*FlankSize)
lines(MedTD_Full, col = "red")
lines(MedTD_Fragm, col = "blue")
points(MedTD_Fragm, col = "blue")


# Plot mean Tajima's D against L1 width
L1Width <- width(L1GRanges[idxUniqueL1Neighbor][unique(OL_L1TD@from)])
InsLengthClass <- cut(L1Width, breaks = seq(0, 7000, 1000))
TDAggregated_mean <- aggregate(data.frame(TD = TDMat[MidPt, ], L1Width = L1Width),
                               by = list(InsLengthClass), 
                               FUN = function(x) mean(x, na.rm = T))

# Check whether there is a relationship between L1 width and Tajima's D among fragments
blnFull <- width(L1GRanges[idxUniqueL1Neighbor][unique(OL_L1TD@from)]) >= 5900
LM1 <- lm(TDMat[MidPt, ] ~ 
          width(L1GRanges[idxUniqueL1Neighbor][unique(OL_L1TD@from)]) + blnFull)
summary(LM1)
LM2 <- lm(TDMat[MidPt, ] ~ 
           width(L1GRanges[idxUniqueL1Neighbor][unique(OL_L1TD@from)]))
summary(LM2)
LM3 <- lm(TDMat[MidPt, !blnFull] ~ 
            width(L1GRanges[idxUniqueL1Neighbor][unique(OL_L1TD@from)])[!blnFull])
summary(LM3)

# Test for different Tajima's D between full-length and fragment at the midpoint
t.test(TDMat[MidPt,idxFull], TDMat[MidPt,idxFragm])

# Trend in mean Tajima's D
LMTrend_Fragm <- lm(MeanTD_Fragm ~ c(1:length(MeanTD_Fragm)))
summary(LMTrend_Fragm)
LMTrend_Full <- lm(MeanTD_Full[-MidPt] ~ c(1:length(MeanTD_Full))[-MidPt])
summary(LMTrend_Full)

# Trend in median Tajima's D
LMMedTrend_Fragm <- lm(MedTD_Fragm[-c(MidPt, NrInt)] ~ 
                         c(1:length(MedTD_Fragm))[-c(MidPt, NrInt)])
summary(LMMedTrend_Fragm)
LMMedTrend_Full <- lm(MedTD_Full[-c(MidPt, NrInt)] ~ 
                        c(1:length(MedTD_Full))[-c(MidPt, NrInt)])
summary(LMMedTrend_Full)

acf(MedTD_Fragm)
acf(MedTD_Full)

# Plot median Tajima's D against absolute difference from center
plot(LMTrend_Full$residuals[(MidPt - 1):1], 
     LMTrend_Full$residuals[(MidPt + 1):(length(MedTD_Full) - 1)])
plot(LMMedTrend_Full$residuals[(MidPt - 1):1], 
     LMMedTrend_Full$residuals[(MidPt + 1):(length(MedTD_Full) - 1)])
cor.test(LMMedTrend_Full$residuals[(MidPt - 1):1], 
     LMMedTrend_Full$residuals[(MidPt + 1):(length(MedTD_Full) - 1)])
cor.test(LMMedTrend_Fragm$residuals[(MidPt - 1):1], 
         LMMedTrend_Fragm$residuals[(MidPt + 1):(length(MedTD_Full) - 1)])
cor.test(MedTD_Full[-MidPt], MedTD_Fragm[-MidPt])
cor.test(LMMedTrend_Full$residuals, LMMedTrend_Fragm$residuals)
plot(LMMedTrend_Full$residuals, LMMedTrend_Fragm$residuals)

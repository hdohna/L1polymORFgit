# The following script compares simulated and estimated L1 and adjusts
# the estimated L1 width to account for detection inaccuracoes

# Load packages
library(GenomicRanges)
library(MASS)
library(fields)
library(raster)
library(lattice)

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')
par(mfrow = c(1, 1))
##########################################
#                                        #
#     Define functions                   #
#                                        #
##########################################

# Function to add columns to L1 detection dataframe
AddCols2L1Detect <- function(L1DetectDF, L1WBreaks = seq(0, 6500, 500)){
  L1DetectDF$L1StartTrue    <- as.numeric(as.character(L1DetectDF$L1StartTrue))
  L1DetectDF$L1EndTrue      <- as.numeric(as.character(L1DetectDF$L1EndTrue))
  L1DetectDF$blnPass        <- L1DetectDF$EstFilter == "PASS"
  L1DetectDF$blnDetectPass       <- L1DetectDF$blnDetect & L1DetectDF$blnPass
  L1DetectDF$L1WidthFromStartEnd <- L1DetectDF$L1EndEst - L1DetectDF$L1StartEst
  L1DetectDF$L1WidthAbsDiff      <- abs(L1DetectDF$L1widthTrue - 
                                        L1DetectDF$L1widthEst)
  L1DetectDF$L1WidthStEAbsDiff   <- abs(L1DetectDF$L1widthTrue - 
                                        L1DetectDF$L1WidthFromStartEnd)
  L1DetectDF$L1StartAbsDiff      <- abs(L1DetectDF$L1StartTrue - 
                                          L1DetectDF$L1StartEst)
  L1DetectDF$L1EndAbsDiff   <- abs(L1DetectDF$L1EndTrue - L1DetectDF$L1EndEst)
  L1DetectDF$L1PosAbsDiff   <- abs(L1DetectDF$PosTrue - L1DetectDF$PosEst)
  L1DetectDF$blnFullTrue    <- L1DetectDF$L1widthTrue >= 6000
  L1DetectDF$blnFullEst     <- L1DetectDF$L1widthEst >= 6000
  L1DetectDF$L1widthTrueBin <- cut(L1DetectDF$L1widthTrue, breaks = L1WBreaks)
  L1DetectDF$L1widthEstBin <- cut(L1DetectDF$L1widthEst, breaks = L1WBreaks)
  L1DetectDF
}

# Function to aggregate L1Detect
AggL1Detect <- function(L1DetectDF, L1PosRange = 30){
  
  # Create an ID per L1 insertion (close insertions are considered one)
  L1DetectDF$L1ID  <- CollapseClosePos_idx(DF = L1DetectDF, 
                                           ChromCol = "Chrom", 
                                           PosCol = "PosEst", 
                                           OLRange = L1PosRange, 
                                           blnPairwise = F)
  
  # Aggregate by L1 ID
  L1DetectAgg <- AggDataFrame(L1DetectDF, 
                              GroupCol = "L1ID", 
                              MeanCols = c("L1widthEst", "L1StartEst", "L1EndEst"),
                              SumCols = c("L1GenoEst", "L1GenoTrue"),
                              MedCols = c("L1widthEst", "L1StartEst", "L1EndEst"),
                              MaxCols = "L1EndEst",
                              MinCols = "L1StartEst",
                              LengthCols = "L1widthEst",
                              Addcols = c("Chrom", "PosTrue", "L1widthTrue", 
                                          "L1StartTrue", "L1EndTrue"))
  
  # Add additional columns
  L1DetectAgg$L1Width_minmax     <- L1DetectAgg$L1EndEst_max - 
    L1DetectAgg$L1StartEst_min
  L1DetectAgg$L1Width_mix        <- L1DetectAgg$L1Width_minmax
  L1DetectAgg$L1Width_mix[L1DetectAgg$L1Width_minmax < 1000] <- 
    L1DetectAgg$L1widthEst_med[L1DetectAgg$L1Width_minmax < 1000]
  L1DetectAgg$L1WidthAbsDiff_med    <- abs(L1DetectAgg$L1widthTrue - 
                                             L1DetectAgg$L1widthEst_med)
  L1DetectAgg$L1WidthAbsDiff_minmax <- abs(L1DetectAgg$L1widthTrue - 
                                             L1DetectAgg$L1Width_minmax)
  L1DetectAgg$L1WidthAbsDiff_mix <- abs(L1DetectAgg$L1widthTrue - 
                                          L1DetectAgg$L1Width_mix)
  L1DetectAgg
  
  }  

##########################################
#                                        #
#            Load data                   #
#                                        #
##########################################

# Load data (created in file CompareSimulatedEstimated L1and add columns
load("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1Simulated_MELT.RData")
L1Detect       <- AddCols2L1Detect(L1Detect)
L1Detect_Group <- AddCols2L1Detect(L1Detect_Group)

# Add numeric columns for L1 start and end
L1_1000G$L1StartNum <- as.numeric(as.character(L1_1000G$L1Start))
L1_1000G$L1EndNum   <- as.numeric(as.character(L1_1000G$L1End))

# Add columns for estimated L1 length to L1DetectAgg_Group
L1DetectAgg_Group$L1widthEst <- sapply(L1DetectAgg_Group$INFO, GetFromVcfINFO_SVLength)
L1StartEnd  <- t(sapply(VcfFile$INFO, GetFromVcfINFO_MELT_L1StartEnd))
L1DetectAgg_Group$L1StartEst <- L1StartEnd[,1]
L1DetectAgg_Group$L1EndEst   <- L1StartEnd[,2]
L1DetectAgg_Group$L1WidthAbsDiff  <- abs(L1DetectAgg_Group$L1widthTrue - 
                                           L1DetectAgg_Group$L1widthEst)

# Add frequency columns
L1DetectAgg_Group$TrueFreq <- rowSums(L1DetectAgg_Group[,as.character(SampleIDs)], na.rm = T)
idxEstGenoCol <- grep("hg19_", colnames(L1DetectAgg_Group))
for (i in idxEstGenoCol){
  L1DetectAgg_Group[,i] <- sapply(L1DetectAgg_Group[,i], GetFromVcfGeno_GenoNum)
}
L1DetectAgg_Group$EstFreq  <- rowSums(L1DetectAgg_Group[, idxEstGenoCol], na.rm = T)

##########################################
#                                        #
#     Plot joint frequencies           #
#                                        #
##########################################

L1widthBreaks <- seq(0, 6500, 500)

# Cut true and estimated L1 width into bins
SimL1widthTrueCut <- cut(L1DetectAgg_Group$L1widthTrue, breaks = L1widthBreaks)
SimL1widthEstCut  <- cut(L1DetectAgg_Group$L1widthEst, breaks  = L1widthBreaks)

# Get combinations of L1 width classes among estimated 
L1combos <- table(SimL1widthTrueCut, SimL1widthEstCut)
levelplot(L1combos,  scales = list(x = list(rot = 90)), xlab = "True LINE-1 width", 
          ylab = "Estimated LINE-1 width")
CreateDisplayPdf('D:/L1ManuscriptFigures/TrueVsEstL1Length.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7, width = 7)

# par(mai = c(2, 2, 0.1, 2))
# image(L1combos, xaxt = "n", yaxt = "n", col = topo.colors(50))
# axis(1, at = seq(0, 1, 1/(length(L1widthBreaks) - 2)), labels = row.names(L1combos),
#      las = 2)
# axis(2, at = seq(0, 1, 1/(length(L1widthBreaks) - 2)), labels = row.names(L1combos),
#      las = 2)
# mtext(text = "Estimated L1 width", side = 1, las = 1, line = 8)
# mtext(text = "Estimated L1 width", side = 1, las = 1, line = 8)
# legend("topright")

##########################################
#                                        #
#     Analyze L1 width estimation        #
#     as function of frequency           #
#                                        #
##########################################

# Fit logistic regression for probability of estimated length within 50 bp of
# true length against estimated frequency
LogReg_L1widthEstvsFreq <- glm(L1WidthAbsDiff < 50 ~ EstFreq, data = L1DetectAgg_Group,
           family = binomial)
summary(LogReg_L1widthEstvsFreq)

# Plot fitted model
FreqVals  <- 1:200
TrueL1exp <- exp(LogReg_L1widthEstvsFreq$coefficients[1] + 
                   LogReg_L1widthEstvsFreq$coefficients[2] * FreqVals)
TrueL1P <- TrueL1exp / (1 + TrueL1exp)
plot(FreqVals, TrueL1P, type = "l",
     xlab = "Sample frequency", ylab = "P(true length estimated)")

##########################################
#                                        #
#     Estimate probability of            #
#  estimated length given true length    #
#                                        #
##########################################

############
# For L1DetectAgg_Group:
###########
L1DetectDF <- L1DetectAgg_Group
L1DetectDF <- L1DetectDF[which(L1DetectDF$L1widthTrue > 0 &
                               L1DetectDF$L1widthEst > 0), ]

plot(L1DetectDF$L1widthTrue, L1DetectDF$L1widthEst)

# Estimate joint surface for true and estimated L1 width
PTrueEstKD <- kde2d(x = L1DetectDF$L1widthTrue, 
                    y = L1DetectDF$L1widthEst,
                    h = 0.1 *bandwidth.nrd(L1DetectDF$L1widthTrue),
      n = 100,
      lims = c(1, 6018, 1, 6019))
image(PTrueEstKD, col = topo.colors(10))
points(L1DetectDF$L1widthTrue, L1DetectDF$L1widthEst, col = "red")
bandwidth.nrd(L1DetectDF$L1widthTrue)

# Calculate estimated L1 length conditional on true length (PEstGivenTrue)
PTrueEst <- PTrueEstKD$z
diag(PTrueEst) <- 0
PEstGivenTrue <- PTrueEst / rowSums(PTrueEst)
image(PEstGivenTrue)

# Create index to map L1widthTrue to the bins of PTrueEstKD
idxClosestX <- sapply(1:nrow(L1DetectDF), function(i) {
  which.min(abs(L1DetectDF$L1widthTrue[i] - PTrueEstKD$x))
})
idxClosestXCount <- table(idxClosestX)
PL1TrueBin <- rep(0, 100)
PL1TrueBin[as.numeric(names(idxClosestXCount))] <- idxClosestXCount

idxClosestY <- sapply(1:nrow(L1DetectDF), function(i) {
  which.min(abs(L1DetectDF$L1widthEst[i] - PTrueEstKD$y))
})
idxClosestYCount <- table(idxClosestY)
PL1EstBin <- rep(0, 100)
PL1EstBin[as.numeric(names(idxClosestYCount))] <- idxClosestYCount
plot(PL1EstBin, type = "l", ylim = c(0, 130))
lines(PL1TrueBin, col = "red")

# # Create matrix of joint probability of true and estimated L1 length
# # (row = true length, col = estimated length)
# PTrueAndEst <- matrix(0, nrow = max(L1DetectDF$L1widthTrue, na.rm = T),
#                       ncol = max(L1DetectDF$L1widthEst, na.rm = T))
# for (i in 1:nrow(L1DetectDF)){
#   idxR <- L1DetectDF$L1widthTrue[i]
#   idxC <- L1DetectDF$L1widthEst[i]
#   PTrueAndEst[idxR, idxC] <- PTrueAndEst[idxR, idxC] + 1
# }
# #image(PTrueAndEst)
# PTrueAndEst <- PTrueAndEst / sum(PTrueAndEst)
# 
# # Create matrix of estimated L1 length conditional on true L1 length
# PEstGivenTrue <- PTrueAndEst / rowSums(PTrueAndEst)




##########################################
#                                        #
#    Calculate summary quantities        #
#                                        #
##########################################

# Compare true and estimated genotype
table(L1Detect$L1GenoTrue, L1Detect$L1GenoEst)
table(L1Detect_Group$L1GenoTrue, L1Detect_Group$L1GenoEst)

# Get proportion of correct L1 width, start, and end estimates
blnGen0 <- L1Detect$L1GenoEst == 0
mean(L1Detect$L1WidthAbsDiff[!blnGen0]    < 100, na.rm = T)
mean(L1Detect$L1WidthStEAbsDiff[!blnGen0] < 100, na.rm = T)
mean(L1Detect$L1StartAbsDiff[!blnGen0]    < 100, na.rm = T)
mean(L1Detect$L1EndAbsDiff[!blnGen0]      < 100, na.rm = T)
mean(L1Detect$L1StartAbsDiff[!blnGen0]    < 100 &
     L1Detect$L1EndAbsDiff[!blnGen0]      < 100, na.rm = T)
mean(L1Detect$L1StartAbsDiff[!blnGen0]    < 100, na.rm = T) *
mean(L1Detect$L1EndAbsDiff[!blnGen0]      < 100, na.rm = T)

blnGen0_G <- L1Detect_Group$L1GenoEst == 0
mean(L1Detect_Group$L1WidthAbsDiff[!blnGen0_G]    < 100, na.rm = T)
mean(L1Detect_Group$L1WidthStEAbsDiff[!blnGen0_G] < 100, na.rm = T)
mean(L1Detect_Group$L1StartAbsDiff[!blnGen0_G]    < 100, na.rm = T)
mean(L1Detect_Group$L1EndAbsDiff[!blnGen0_G]      < 100, na.rm = T)
mean(L1Detect_Group$L1StartAbsDiff[!blnGen0_G]    < 100 &
       L1Detect_Group$L1EndAbsDiff[!blnGen0_G]    < 100, na.rm = T)
mean(L1Detect_Group$L1StartAbsDiff[!blnGen0_G]    < 100, na.rm = T) *
  mean(L1Detect_Group$L1EndAbsDiff[!blnGen0_G]    < 100, na.rm = T)

mean(L1DetectAgg$L1WidthAbsDiff < 100, na.rm = T)
mean(L1DetectAgg$L1WidthAbsDiff_med    < 100)
mean(L1DetectAgg$L1WidthAbsDiff_minmax < 100)
mean(L1DetectAgg$L1WidthAbsDiff_mix < 100)

# Proportion of estimated full-length that are full-length
blnTrueFull <- L1Detect$L1widthTrue  >= 6000
blnEstFull  <- L1Detect$L1widthEst   >= 6000
blnEndFull  <- L1Detect$L1EndEst     >= 5500
blnStartFull  <- L1Detect$L1StartEst <= 500
# (i) all (single)
sum(blnTrueFull & blnEstFull & !blnGen0, na.rm = T) / 
  sum(blnEstFull & !blnGen0, na.rm = T)
# (ii) among L1 with 3' end covered
sum(blnTrueFull & blnEstFull & blnEndFull & !blnGen0, na.rm = T) / 
  sum(blnEstFull & blnEndFull & !blnGen0, na.rm = T)
# (iii) among L1 with 5' and 3' end covered
sum(blnTrueFull & blnEstFull & blnEndFull & blnStartFull & !blnGen0, na.rm = T) / 
  sum(blnEstFull & blnEndFull & blnStartFull & !blnGen0, na.rm = T)
# (iv) among L1 with 5' end not covered  and 3' end covered
sum(blnTrueFull & blnEstFull & blnEndFull & (!blnStartFull) & !blnGen0, na.rm = T) / 
  sum(blnEstFull & blnEndFull & (!blnStartFull) & !blnGen0, na.rm = T)

plot(L1Detect$L1EndEst[blnEstFull & !blnGen0], 
     L1Detect$L1widthTrue[blnEstFull & !blnGen0])
plot(L1Detect$L1StartEst[blnEstFull], L1Detect$L1widthTrue[blnEstFull])

# Proportion of estimated full-length that are full-length
# (i) all (single)
sum(L1Detect$L1widthTrue  >= 6000 & L1Detect$L1widthEst >= 6000, na.rm = T) /
  sum(L1Detect$L1widthEst >= 6000, na.rm = T)
# (ii) PASS only (single)
sum(L1Detect$L1widthTrue  >= 6000 & L1Detect$L1widthEst >= 6000 
    & L1Detect$blnPass, na.rm = T) /
  sum(L1Detect$L1widthEst >= 6000 & L1Detect$blnPass, na.rm = T)
# (iii) all (group)
sum(L1Detect_Group$L1widthTrue >= 6000 & L1Detect_Group$L1widthEst >= 6000 &
      !blnGen0_G, na.rm = T) /
  sum(L1Detect_Group$L1widthEst >= 6000 & !blnGen0_G, na.rm = T)
# (iv) PASS only (group)
sum(L1Detect_Group$L1widthTrue >= 6000 & 
      L1Detect_Group$L1widthEst >= 6000 
    & L1Detect_Group$blnPass & !blnGen0_G, na.rm = T) /
  sum(L1Detect_Group$L1widthEst >= 6000 & L1Detect_Group$blnPass &
        !blnGen0_G, na.rm = T)

# Proportion of short L1s that are labelled as full-length
sum(L1Detect$L1widthTrue <= 1000 & L1Detect$L1widthEst >= 6000, na.rm = T) /
  sum(L1Detect$L1widthTrue <= 1000, na.rm = T)

# Get info of L1s that are estimated to be full-length but are not
L1Detect$L1EstINFO[which(L1Detect$L1widthTrue <= 1000 & 
                           L1Detect$L1widthEst >= 6000)]
L1FalseFull <- L1Detect[which(L1Detect$L1widthTrue <= 1000 & 
                                L1Detect$L1widthEst >= 6000),
                c("L1StartTrue", "L1EndTrue", "L1widthTrue", "L1widthEst", "L1EstINFO")]

# Sensitivity and specificity
Sensitivity_all <- mean(L1Detect$blnDetectPass, na.rm = T)
blnL1Estimated  <- !is.na(L1Detect$L1GenoEst)
Specificity_all <- sum(blnL1Estimated & L1Detect$blnDetectPass, na.rm = T) / 
  sum(blnL1Estimated)
cat("Sensitivity:", Sensitivity_all, "\n")
cat("Specificity:", Specificity_all, "\n")

# Sensitivity and specificity in group
Sensitivity_all <- mean(L1Detect_Group$blnDetectPass, na.rm = T)
blnL1Estimated  <- !is.na(L1Detect_Group$L1GenoEst)
Specificity_all <- sum(blnL1Estimated & L1Detect_Group$blnDetectPass, na.rm = T) / 
  sum(blnL1Estimated)
cat("Group sensitivity:", Sensitivity_all, "\n")
cat("Group specificity:", Specificity_all, "\n")

# Proportion that over -and underestimates L1 length
mean(L1Detect$L1widthTrue < L1Detect$L1widthEst - 100, na.rm = T)
mean(L1Detect$L1widthTrue > L1Detect$L1widthEst + 100 , na.rm = T)

##########################################
#                                        #
#         Analyze frequencies            #
#                                        #
##########################################

# Get unique IDs
UniqueIDs   <- unique(as.character(L1Detect$SampleID))
UniqueIDs_G <- unique(as.character(L1Detect_Group$SampleID))
length(UniqueIDs)
length(UniqueIDs_G)

# # Get frequency from 1000 genome data
# L1Freq1KG <- rowSums(L1_1000G[,UniqueIDs])
# 
# # Subset to get
# L1DetectAggTrue <- L1DetectAgg[!is.na(L1DetectAgg$PosTrue), ]
# 
# # Match 1000 genome entries to L1 detect
# L1DetectAgg_GR <- makeGRangesFromDataFrame(L1DetectAggTrue,
#                                             seqnames.field = "Chrom",
#                                             start.field = "PosTrue",
#                                             end.field = "PosTrue")
# OL_1000G_Detect <- findOverlaps(L1_1000G_GR_hg19, L1DetectAgg_GR)
# 
# # Append true frequency 
# L1DetectAggTrue$FreqTrue <- NA
# L1DetectAggTrue$FreqTrue[OL_1000G_Detect@to] <- 
#   L1Freq1KG[OL_1000G_Detect@from]
# 
# plot(L1DetectAggTrue$FreqTrue, L1DetectAggTrue$L1GenoEst_sum)
# lines(c(0, 60), c(0, 60), col = "red")

plot(L1DetectAgg$L1GenoTrue_sum, L1DetectAgg$L1GenoEst_sum)
lines(c(0, 60), c(0, 60), col = "red")

# Compare true and estimated mean frequency per length class
L1WBreaks <- seq(0, 7000, 1000)
L1DetectAgg$L1widthClassTrue <- cut(L1DetectAgg$L1widthTrue, breaks = L1WBreaks)
L1DetectAgg$L1widthClassEst  <- cut(L1DetectAgg$L1widthEst_med, breaks = L1WBreaks)
L1DetectAggPerL1widthTrue <- AggDataFrame(L1DetectAgg, 
                                                GroupCol = "L1widthClassTrue", 
                                                MeanCols = c("L1widthTrue", "L1GenoTrue_sum"))
L1DetectAggPerL1widthEst <- AggDataFrame(L1DetectAgg, 
                                               GroupCol = "L1widthClassEst", 
                                               MeanCols = c("L1widthEst_med", "L1GenoEst_sum"))
plot(L1DetectAggPerL1widthEst$L1widthEst_med_mean, 
     L1DetectAggPerL1widthEst$L1GenoEst_sum_mean, col = "red")
points(L1DetectAggPerL1widthTrue$L1widthTrue_mean, 
       L1DetectAggPerL1widthTrue$L1GenoTrue_sum_mean, col = "blue")
plot(L1DetectAggPerL1widthEst$L1widthEst_med_mean, 
     L1DetectAggPerL1widthEst$L1GenoEst_sum_mean - 
      L1DetectAggPerL1widthTrue$L1GenoTrue_sum_mean)
lines(c(0, 10^4), c(0, 0), lty =2) 
 

L1DetectAgg$blnFull <- L1DetectAgg$L1widthTrue >= 6000
L1DetectAgg$FreqDiff <- L1DetectAgg$RelFreq - L1DetectAgg$FreqTrue
L1DetectAgg$RelFreqDiff <- (L1DetectAgg$RelFreq - L1DetectAgg$FreqTrue)/
  L1DetectAgg$FreqTrue

LM_BiasL1Width <- lm(FreqDiff ~ L1widthTrue + blnFull, data = L1DetectAgg)
summary(LM_BiasL1Width)

LM_RelBiasL1Width <- lm(RelFreqDiff ~ L1widthTrue, data = L1DetectAgg)
summary(LM_BiasL1Width)

# Check whether among true full-length L1s the estimated length increases
# with frequency
LM_L1widthEst_med_VsFreq <- lm(L1widthEst_med ~ L1GenoEst_sum,
                               data = L1DetectAgg[L1DetectAgg$L1widthTrue >= 6000, ])
summary(LM_L1widthEst_med_VsFreq)
LM_L1widthEst_minmax_VsFreq <- lm(L1Width_minmax ~ L1GenoEst_sum,
                                  data = L1DetectAgg[L1DetectAgg$L1widthTrue >= 6000, ])
summary(LM_L1widthEst_minmax_VsFreq)
LM_L1widthEst_mix_VsFreq <- lm(L1Width_mix ~ L1GenoEst_sum,
                               data = L1DetectAgg[L1DetectAgg$L1widthTrue >= 6000, ])
summary(LM_L1widthEst_mix_VsFreq)

# Check whether among true full-length L1s the estimated length increases
# with frequency
L1DetectAgg_Group$TrueFreq <- rowSums(L1DetectAgg_Group[,as.character(SampleIDs)], na.rm = T)
idxEstGenoCol <- grep("hg19_", colnames(L1DetectAgg_Group))
for (i in idxEstGenoCol){
  L1DetectAgg_Group[,i] <- sapply(L1DetectAgg_Group[,i], GetFromVcfGeno_GenoNum)
}
L1DetectAgg_Group$EstFreq  <- rowSums(L1DetectAgg_Group[, idxEstGenoCol], na.rm = T)
LM_L1widthEst_VsFreq <- lm(L1widthEst ~ TrueFreq,
                           data = L1DetectAgg_Group[L1DetectAgg_Group$L1widthTrue >= 6000, ])
summary(LM_L1widthEst_VsFreq)
plot(L1DetectAgg_Group$TrueFreq, 0.8*L1DetectAgg_Group$EstFreq)
lines(c(0, 100), c(0, 100))
L1DetectAgg_Group[L1DetectAgg_Group$EstFreq > 80, ]
L1DetectAgg_Group$EstFreqAdj <- 0.8*L1DetectAgg_Group$EstFreq

# Compare true and estimated mean frequency per length class
L1WBreaks <- seq(0, 7000, 1000)
L1DetectAgg_Group$L1widthClassTrue <- cut(L1DetectAgg_Group$L1widthTrue, breaks = L1WBreaks)
L1DetectAgg_Group$L1widthClassEst  <- cut(L1DetectAgg_Group$L1widthEst, breaks = L1WBreaks)
L1DetectAgg_GroupPerL1widthTrue <- AggDataFrame(L1DetectAgg_Group[L1DetectAgg_Group$EstFreq < 86,], 
                            GroupCol = "L1widthClassTrue", 
                            MeanCols = c("L1widthTrue", "TrueFreq"))
L1DetectAgg_GroupPerL1widthEst <- AggDataFrame(L1DetectAgg_Group[L1DetectAgg_Group$EstFreq < 86,], 
                                             GroupCol = "L1widthClassEst", 
                                             MeanCols = c("L1widthEst", "EstFreq", "EstFreqAdj"))
plot(L1DetectAgg_GroupPerL1widthEst$L1widthEst_mean, 
     L1DetectAgg_GroupPerL1widthEst$EstFreq_mean, col = "red")
points(L1DetectAgg_GroupPerL1widthTrue$L1widthTrue_mean, 
       L1DetectAgg_GroupPerL1widthTrue$TrueFreq_mean, col = "blue")
plot(L1DetectAgg_GroupPerL1widthEst$L1widthEst_mean, 
     L1DetectAgg_GroupPerL1widthEst$EstFreq_mean - 
       L1DetectAgg_GroupPerL1widthTrue$TrueFreq_mean)
plot(L1DetectAgg_GroupPerL1widthEst$L1widthEst_mean, 
     L1DetectAgg_GroupPerL1widthEst$EstFreqAdj_mean - 
       L1DetectAgg_GroupPerL1widthTrue$TrueFreq_mean)


##########################################
#                                        #
#        Logistic regression             #
#          for detection                 #
#                                        #
##########################################

# Detection as function of predictors
LogReg_DetectL1widthTrue <- glm(blnDetectPass ~ L1widthTrue + blnFullTrue, 
                                data = L1Detect[L1Detect$L1GenoTrue == 1,],
                               family = binomial)
summary(LogReg_DetectL1widthTrue)

LogReg_DetectL1StETrue <- glm(blnDetect ~ L1StartTrue + L1widthTrue +
                                blnFullTrue, data = L1Detect,
                                family = binomial)
summary(LogReg_DetectL1StETrue)
LogReg_DetectL1widthTrueSum <- summary(LogReg_DetectL1widthTrue)
LogReg_DetectL1widthTrue$coefficients
LogReg_DetectL1widthEst <- glm(blnDetect ~ L1StartEst, 
                               data = L1Detect,
                               family = binomial)
summary(LogReg_DetectL1widthEst)


# # Plot estimated start vs estimated start from discorant read pairs
# plot  (L1Detect$L1StartEst, L1Detect$L1StartEst_DiscReads,
#        xlab = "Estimated L1 start", ylab = "L1 start estimated from disc. reads")
# lines (c(0, 6000), c(0, 6000), col = "red", lwd = 2)
# 
# # Plot estimated start vs estimated start from discorant read pairs
# plot  (L1Detect$L1StartTrueNum, L1Detect$L1StartEst_DiscReads,
#        xlab = "True L1 start", ylab = "L1 start estimated from disc. reads")
# lines (c(0, 6000), c(0, 6000), col = "red", lwd = 2)

plot  (L1Detect_Group$L1widthTrue, L1Detect_Group$L1widthEst)
hist  (L1Detect$L1widthTrue)
hist  (L1Detect$L1widthEst)
table (L1Detect$L1widthEst)
CurrentChrom = "chr1"
ChrT <- CreateInsertTxt(1:5)

L1Detect[which(L1Detect$SampleID == "HG00107")[1:10], ]
L1_1000G[1:5,  c("CHROM", "POS", "L1Start", "L1End")]

L1Detect$L1StartEst[1:5]
as.numeric(L1Detect$L1StartEst)[1:5]
L1Detect$L1StartTrue[1:5]
as.numeric(as.character(L1Detect$L1StartTrue))[1:5]

plot(as.numeric(as.character(L1Detect$L1StartTrue)), as.numeric(L1Detect$L1StartEst))

par(mfrow = c(1, 2), oma = c(3, 4, 1, 1),
    mar = )
plot(1:10, 1:10, ylab = "", xlab = "X label")
plot(1:10, 1:10, ylab = "", xlab = "X label")
mtext(side = 2, "Y label", outer = T)


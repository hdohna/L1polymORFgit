# The following script compares simulated and estimated L1.
# It uses objects created in script 'CompareSimulatedEstimatedL1.R'

# Load packages
library(GenomicRanges)

##########################################
#                                        #
#     Define functions                   #
#                                        #
##########################################

# Function to add columns to L1 detection dataframe
AddCols2L1Detect <- function(L1DetectDF){
  L1DetectDF$L1StartTrue    <- as.numeric(as.character(L1DetectDF$L1StartTrue))
  L1DetectDF$L1EndTrue      <- as.numeric(as.character(L1DetectDF$L1EndTrue))
  L1DetectDF$blnPass        <- L1DetectDF$EstFilter == "PASS"
  L1DetectDF$L1WidthFromStartEnd <- L1DetectDF$L1EndEst - L1DetectDF$L1StartEst
  L1DetectDF$L1WidthAbsDiff <- abs(L1DetectDF$L1widthTrue - L1DetectDF$L1widthEst)
  L1DetectDF$L1WidthStEAbsDiff <- abs(L1DetectDF$L1widthTrue - L1DetectDF$L1WidthFromStartEnd)
  L1DetectDF$L1StartAbsDiff <- abs(L1DetectDF$L1StartTrue - L1DetectDF$L1StartEst)
  L1DetectDF$L1EndAbsDiff   <- abs(L1DetectDF$L1EndTrue - L1DetectDF$L1EndEst)
  L1DetectDF$L1PosAbsDiff   <- abs(L1DetectDF$PosTrue - L1DetectDF$PosEst)
  L1DetectDF$blnFullTrue    <- L1DetectDF$L1widthTrue >= 6000
  L1DetectDF$blnFullEst     <- L1DetectDF$L1widthEst >= 6000
  L1DetectDF
}


##########################################
#                                        #
#            Load data                   #
#                                        #
##########################################

# Load data and add columns
load("D:/L1polymORF/Data/L1Simulated_MELT.RData")
L1Detect       <- AddCols2L1Detect(L1Detect)
L1Detect_Group <- AddCols2L1Detect(L1Detect_Group)

which(substr(L1Detect$Chrom, 1, 3) == "chr")
L1Detect[163,]

# Add numeric columns for L1 start and end
L1_1000G$L1StartNum <- as.numeric(as.character(L1_1000G$L1Start))
L1_1000G$L1EndNum   <- as.numeric(as.character(L1_1000G$L1End))

##########################################
#                                        #
#          Aggregate detection           #
#            per L1 insertion            #
#                                        #
##########################################

# Specifiy a range within insertions are considered to be the same
L1PosRange = 30

# Subset to retain only entries with estimated L1 position
L1DetectDF <- L1Detect[!is.na(L1Detect$PosEst),]

# Create an ID per L1
L1DetectDF$L1ID <- paste(L1DetectDF$Chrom, L1DetectDF$PosEst)
blnDupl    <- duplicated(L1DetectDF$L1ID)
L1IDUnique <- L1DetectDF$L1ID[!blnDupl]

# Create genomic ranges for unique positions
L1Detect_GR <- makeGRangesFromDataFrame(L1DetectDF[!blnDupl, ],
                                        seqnames.field = "Chrom",
                                        start.field = "PosEst",
                                        end.field = "PosEst")
# Resize genomic ranges for unique positions and find overlaps with
# itself
L1Detect_GR_large <- resize(L1Detect_GR, 2*L1PosRange, fix = "center")
L1Detect_OL <- findOverlaps(L1Detect_GR_large, L1Detect_GR_large)

# Assign to all overlapping ranges the smallest index
L1idx <- sapply(1:length(L1IDUnique), function(x){
           blnOL <- L1Detect_OL@from ==x | L1Detect_OL@to == x
           min(c(L1Detect_OL@from[blnOL], L1Detect_OL@to[blnOL]))
         })
L1idxMatch  <- match(L1DetectDF$L1ID, L1IDUnique)
L1DetectDF$L1ID <- L1IDUnique[L1idx][L1idxMatch]

# Remove "chr" from IDs that start with "chr" [FIGURE OUT WHY]
blnChr <- substr(L1DetectDF$L1ID, 1, 3) == "chr"
any(blnChr)
L1DetectDF$L1ID[blnChr] <- substr(L1DetectDF$L1ID[blnChr], 4, 
                                  nchar(L1DetectDF$L1ID[blnChr]))

# Get L1 frequency, mean and medium start, end and width
L1DetectAgg_Geno <- aggregate(L1DetectDF$L1GenoEst, by = list(L1DetectDF$L1ID), 
                              FUN = sum)
colnames(L1DetectAgg_Geno) <- c("L1ID", "Freq")
L1DetectAgg_StartMean <- aggregate(L1DetectDF$L1StartEst, by = list(L1DetectDF$L1ID), 
                                   FUN = mean)
colnames(L1DetectAgg_StartMean) <- c("L1ID", "StartMean")
L1DetectAgg_StartMed <- aggregate(L1DetectDF$L1StartEst, by = list(L1DetectDF$L1ID), 
                                  FUN = median)
colnames(L1DetectAgg_StartMed) <- c("L1ID", "StartMed")
L1DetectAgg_EndMean <- aggregate(L1DetectDF$L1EndEst, by = list(L1DetectDF$L1ID), 
                                 FUN = mean)
colnames(L1DetectAgg_EndMean) <- c("L1ID", "EndMean")
L1DetectAgg_EndMed <- aggregate(L1DetectDF$L1EndEst, by = list(L1DetectDF$L1ID), 
                                FUN = median)
colnames(L1DetectAgg_EndMed) <- c("L1ID", "EndMed")
L1DetectAgg_EndMax <- aggregate(L1DetectDF$L1EndEst, by = list(L1DetectDF$L1ID), 
                                FUN = max)
colnames(L1DetectAgg_EndMax) <- c("L1ID", "EndMax")
L1DetectAgg_widthMean <- aggregate(L1DetectDF$L1widthEst, by = list(L1DetectDF$L1ID), 
                                   FUN = mean)
colnames(L1DetectAgg_widthMean) <- c("L1ID", "widthMean")
L1DetectAgg_widthMed <- aggregate(L1DetectDF$L1widthEst, by = list(L1DetectDF$L1ID), 
                                  FUN = median)
colnames(L1DetectAgg_widthMed) <- c("L1ID", "widthMed")

L1DetectAgg <- merge(L1DetectAgg_Geno, L1DetectAgg_StartMean)
L1DetectAgg <- merge(L1DetectAgg, L1DetectAgg_StartMed)
L1DetectAgg <- merge(L1DetectAgg, L1DetectAgg_EndMean)
L1DetectAgg <- merge(L1DetectAgg, L1DetectAgg_EndMed)
L1DetectAgg <- merge(L1DetectAgg, L1DetectAgg_widthMean)
L1DetectAgg <- merge(L1DetectAgg, L1DetectAgg_widthMed)
L1DetectAgg <- merge(L1DetectAgg, L1DetectAgg_EndMax)
L1DetectAgg$Chrom  <- sapply(L1DetectAgg$L1ID, function(x) strsplit(x, " ")[[1]][1])
L1DetectAgg$Pos    <- sapply(L1DetectAgg$L1ID, function(x) as.numeric(strsplit(x, " ")[[1]][2]))
L1IDmatch          <- match(L1DetectAgg$L1ID, L1DetectDF$L1ID)
L1DetectAgg$L1widthEst <- L1DetectDF$L1widthEst[L1IDmatch]
L1DetectAgg$L1widthTrue <- L1DetectDF$L1widthTrue[L1IDmatch]
L1DetectAgg$L1StartTrue <- L1DetectDF$L1StartTrue[L1IDmatch]
L1DetectAgg$L1EndTrue <- L1DetectDF$L1EndTrue[L1IDmatch]

L1DetectAgg$L1WidthAbsDiff <- abs(L1DetectAgg$L1widthTrue - L1DetectAgg$widthMed)

##########################################
#                                        #
#         Overview plots                 #
#                                        #
##########################################

# Plot estimated and true start and end against each other
plot(L1Detect$L1StartTrue, L1Detect$L1StartEst)
plot(L1Detect$L1EndTrue,   L1Detect$L1EndEst)
plot(L1Detect_Group$L1StartTrue, L1Detect_Group$L1StartEst)
plot(L1Detect_Group$L1EndTrue,   L1Detect_Group$L1EndEst)
plot(L1DetectAgg$L1StartTrue, L1DetectAgg$StartMean)
plot(L1DetectAgg$L1StartTrue, L1DetectAgg$StartMed)
plot(L1DetectAgg$L1EndTrue, L1DetectAgg$EndMed)
plot(L1DetectAgg$L1EndTrue, L1DetectAgg$EndMean)
plot(L1DetectAgg$L1EndTrue, L1DetectAgg$EndMax)
plot  (L1DetectAgg$L1widthTrue, L1DetectAgg$widthMean)
plot  (L1DetectAgg$L1widthTrue, L1DetectAgg$widthMed)

plot(L1Detect$L1EndTrue,   L1Detect$L1EndEst)
plot(L1Detect$L1EndEst - L1Detect$L1StartEst, L1Detect$L1widthEst)

# Plot true vs estimated L1 length
plot  (L1Detect$L1widthTrue, L1Detect$L1widthEst)
lines (c(0, 6000), c(0, 6000), col = "red", lwd = 2)
plot  (L1Detect_Group$L1widthTrue, L1Detect_Group$L1widthEst)
lines (c(0, 6000), c(0, 6000), col = "red", lwd = 2)

# Plot histogram in absolute difference between true and estimated insertion
# position
hist(L1Detect$L1PosAbsDiff, breaks = 0:30)
max(L1Detect$L1PosAbsDiff, na.rm = T)
hist(L1Detect_Group$L1PosAbsDiff, breaks = 0:50)
max(L1Detect_Group$L1PosAbsDiff, na.rm = T)

##########################################
#                                        #
#         Calculate quantities           #
#                                        #
##########################################

# Compare true and estimated genotype
table(L1Detect$L1GenoTrue, L1Detect$L1GenoEst)
table(L1Detect_Group$L1GenoTrue, L1Detect_Group$L1GenoEst)

# Get proportion of correct estimates
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
L1Detect$L1EstINFO[which(L1Detect$L1widthTrue <= 1000 & L1Detect$L1widthEst >= 6000)]
L1FalseFull <- L1Detect[which(L1Detect$L1widthTrue <= 1000 & L1Detect$L1widthEst >= 6000),
                        c("L1StartTrue", "L1EndTrue", "L1widthTrue", "L1widthEst", "L1EstINFO")]

# Sensitivity and specificity
Sensitivity_all <- mean(L1Detect$blnDetect * (!blnGen0), na.rm = T)
blnL1Estimated  <- !is.na(L1Detect$L1GenoEst)
Specificity_all <- sum(blnL1Estimated & L1Detect$blnDetect, na.rm = T) / 
  sum(blnL1Estimated)
cat("Sensitivity:", Sensitivity_all, "\n")
cat("Specificity:", Specificity_all, "\n")

# Sensitivity and specificity in group
Sensitivity_all <- mean(L1Detect_Group$blnDetect * (!blnGen0_G), na.rm = T)
blnL1Estimated  <- !is.na(L1Detect_Group$L1GenoEst) & !blnGen0_G
Specificity_all <- sum(blnL1Estimated & L1Detect_Group$blnDetect, na.rm = T) / 
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

# Get frequency from 1000 genome data
L1Freq1KG <- rowSums(L1_1000G[,UniqueIDs])

# Match 1000 genome entries to L1 detect
L1DetectAgg
L1DetectAgg$Chromosome <- paste("chr", L1DetectAgg$Chrom, sep = "")
L1DetectAgg_GR <- makeGRangesFromDataFrame(L1DetectAgg,
                                            seqnames.field = "Chromosome",
                                            start.field = "Pos",
                                            end.field = "Pos")
OL_1000G_Detect <- findOverlaps(L1_1000G_GR_hg19, L1DetectAgg_GR)

# Append true frequency and true insertion length
L1DetectAgg$FreqTrue <- NA
L1DetectAgg$FreqTrue[OL_1000G_Detect@to] <- 
  L1_1000G$Frequency[OL_1000G_Detect@from]

L1DetectAgg$RelFreq <- L1DetectAgg$Freq / 2 / length(UniqueIDs)
plot(L1DetectAgg$FreqTrue, L1DetectAgg$RelFreq)
lines(c(0, 1), c(0, 1))

L1DetectAgg$blnFull <- L1DetectAgg$L1widthTrue >= 6000
L1DetectAgg$FreqDiff <- L1DetectAgg$RelFreq - L1DetectAgg$FreqTrue
L1DetectAgg$RelFreqDiff <- (L1DetectAgg$RelFreq - L1DetectAgg$FreqTrue)/
  L1DetectAgg$FreqTrue

LM_BiasL1Width <- lm(FreqDiff ~ L1widthTrue + blnFull, data = L1DetectAgg)
summary(LM_BiasL1Width)

LM_RelBiasL1Width <- lm(RelFreqDiff ~ L1widthTrue, data = L1DetectAgg)
summary(LM_BiasL1Width)
##########################################
#                                        #
#        Logistic regression             #
#          for detection                 #
#                                        #
##########################################

# Detection as function of predictors
LogReg_DetectL1widthTrue <- glm(blnDetect ~ L1widthTrue + blnFullTrue, 
                                data = L1Detect,
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


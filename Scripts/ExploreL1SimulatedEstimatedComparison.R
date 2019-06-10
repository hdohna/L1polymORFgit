load("D:/L1polymORF/Data/L1Simulated_MELT.RData")

# Add numeric columns for L1 start and end
L1_1000G$L1StartNum <- as.numeric(as.character(L1_1000G$L1Start))
L1_1000G$L1EndNum   <- as.numeric(as.character(L1_1000G$L1End))

# Add additional columns to L1Detect
L1Detect$L1StartTrueNum <- as.numeric(as.character(L1Detect$L1StartTrue))
L1Detect$L1EndTrueNum   <- as.numeric(as.character(L1Detect$L1EndTrue))
L1Detect$blnPass        <- L1Detect$EstFilter == "PASS"


plot  (L1_1000G$InsLength, L1_1000G$L1EndNum  - L1_1000G$L1StartNum)
lines (c(0, 6000), c(0, 6000), col = "red", lwd = 2)
plot  (L1Detect$L1widthTrue, L1Detect$L1widthEst)
lines (c(0, 6000), c(0, 6000), col = "red", lwd = 2)

# Plot estimated start vs estimated start from discorant read pairs
plot  (L1Detect$L1StartEst, L1Detect$L1StartEst_DiscReads,
       xlab = "Estimated L1 start", ylab = "L1 start estimated from disc. reads")
lines (c(0, 6000), c(0, 6000), col = "red", lwd = 2)

# Plot estimated start vs estimated start from discorant read pairs
plot  (L1Detect$L1StartTrueNum, L1Detect$L1StartEst_DiscReads,
       xlab = "True L1 start", ylab = "L1 start estimated from disc. reads")
lines (c(0, 6000), c(0, 6000), col = "red", lwd = 2)

# Proportion of L1s with properly estimated length
L1Length_DiffTrueEst <- abs(L1Detect$L1widthTrue - L1Detect$L1widthEst)
sum(L1Length_DiffTrueEst < 100, na.rm = T) / sum(!is.na(L1Length_DiffTrueEst))

# Proportion of short L1s that are labelled as full-length
sum(L1Detect$L1widthTrue <= 1000 & L1Detect$L1widthEst >= 6000, na.rm = T) /
  sum(L1Detect$L1widthTrue <= 1000, na.rm = T)

# Proportion of estimated full-length that are full-length
# (i) all
sum(L1Detect$L1widthTrue >= 6000 & L1Detect$L1widthEst >= 6000, na.rm = T) /
  sum(L1Detect$L1widthEst >= 6000, na.rm = T)
# (ii) PASS only
sum(L1Detect$L1widthTrue >= 6000 & L1Detect$L1widthEst >= 6000 
    & blnPass, na.rm = T) /
  sum(L1Detect$L1widthEst >= 6000 & blnPass, na.rm = T)


# Logistic regression for detection probability as finction of insertion length
LogReg_DetectL1width <- glm(blnDetect ~ L1widthTrue + blnFull, data = L1Detect,
                            family = binomial)
summary(LogReg_DetectL1width)

# Return some summaries
Sensitivity_all <- mean(L1Detect$blnDetect, na.rm = T)
blnL1Estimated  <- !is.na(L1Detect$L1GenoEst)
Specificity_all <- sum(blnL1Estimated & L1Detect$blnDetect, na.rm = T) / 
  sum(blnL1Estimated)
cat("Sensitivity:", Sensitivity, "\n")
cat("Specificity:", Specificity, "\n")

# Proportion that over -and underestimates L1 length
mean(L1Detect$L1widthTrue < L1Detect$L1widthEst - 100, na.rm = T)
mean(L1Detect$L1widthTrue > L1Detect$L1widthEst + 100 , na.rm = T)


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


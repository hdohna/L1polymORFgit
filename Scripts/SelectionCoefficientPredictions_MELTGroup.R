# The script below estimates selection coefficients of L1 from the 
# 1000 genome data

##########################################
#                                        #
#     Load packages                      #
#                                        #
##########################################

# Source start script
source('D:/OneDrive - American University of Beirut/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(pracma)

##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

# Specify file paths
InputPath <- "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1SelectionResults_MELT_GroupwithSim.RData"

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("\n\nLoading and processing data ...")

# Load previously generated objects
load(InputPath)
cat("done!\n")

###################################################
#                                                 #
#          Pre-process data                       #
#                                                 #
###################################################

# Create a matrix of predictor variables (L1 width and boolean variable for
# for full-length L1)
XMat  <- as.matrix(cbind(1, L1TotData[, c("blnFull", "L1width")]))

# Create a vector of selection coefficients per L1
L1TotData$sVals <- XMat %*% ModelFit_pracma$ML_abc$par
L1TotData$sVals2 <- XMat %*% (2 * ModelFit_pracma$ML_abc$par)

# Create a vector of L1 length classes
L1TotData$InsLengthClass <- cut(L1TotData$L1width, breaks = 
                                  seq(0, 6500, 500))

# Get mean L1 frequency per length
L1WidthAggregated <- AggDataFrame(L1TotData[L1TotData$blnIns,], 
                                  GroupCol = "InsLengthClass", 
                                  MeanCols = c("L1width", "Freq", "blnFull", "sVals",
                                               "sVals2"), 
                                  LengthCols = "Freq",
                                  VarCols = "Freq")
idxL1LengthMatch <- match(L1TotData$InsLengthClass, 
                          L1WidthAggregated$InsLengthClass)                                              

max(1 / L1TotData$L1Freq)
SSize      <- L1TotData$SampleSize[1]


###################################################
#                                                 #
#   Plot frequency vs. insertion length           #
#                                                 #
###################################################

# Get sample size and create a range of s-values
LengthVals <- seq(0, 6400, 200)
Full      <- LengthVals >= 6000
SVals <- ModelFit_pracma$ML_abc$par[1] + ModelFit_pracma$ML_abc$pa[2]*Full +
  ModelFit_pracma$ML_abc$par[3] * LengthVals
SValsTa <- ML_L1WidthFullTa_nonTa$par[1] + ML_L1WidthFullTa_nonTa$pa[2]*Full +
  ML_L1WidthFullTa_nonTa$par[4] * LengthVals
SValsnonTa <-ML_L1WidthFullTa_nonTa$par[1] + ML_L1WidthFullTa_nonTa$pa[3]*Full +
  ML_L1WidthFullTa_nonTa$par[4] * LengthVals
DetectProb <- L1TotData$DetectProb[1]
ML_L1WidthFullTa_nonTa
plot(SVals, SValsTa)
plot(SVals, SValsnonTa)
lines(c(-10, 10), c(-10, 10))
# Calculate expected frequency per L1 width
cat("\nCalculate expected frequency per L1 width ...")
ExpL1Width <- sapply(1:length(SVals), function(i) {
  ExpAlleleFreq_pracma(s = SVals[i], N = PopSize, SampleSize = SSize,
                       DetectProb = DetectProb, blnIns = T, 
                       LogRegCoeff = LogRegL1Ref$coefficients)
})
ExpL1WidthTa <- sapply(1:length(SVals), function(i) {
  ExpAlleleFreq_pracma(s = SValsTa[i], N = PopSize, SampleSize = SSize,
                DetectProb = DetectProb, blnIns = T, 
                LogRegCoeff = LogRegL1Ref$coefficients)
  })
ExpL1Width_nonTa <- ExpL1WidthTa
ExpL1Width_nonTa[Full] <- sapply(which(Full), function(i) {
  ExpAlleleFreq_pracma(s = SValsnonTa[i], N = PopSize, SampleSize = SSize,
                       DetectProb = DetectProb, blnIns = T, 
                       LogRegCoeff = LogRegL1Ref$coefficients)
})
cat("done!\n")

# Plot mean frequency, expected frequency
par(mfrow = c(1, 1), oma = c(0.2,  0.2,  0.2,  0.2), 
     mai = c(1, 1, 0.2, 1),
     cex.lab = 1)
plot(L1WidthAggregated$L1width_mean, 
     L1WidthAggregated$Freq_mean/SSize, xlab = "LINE-1 length [bp]",
     ylab = "Mean LINE-1 frequency", ylim = c(0, 0.03))
AddErrorBars(MidX = L1WidthAggregated$L1width_mean, 
             MidY = L1WidthAggregated$Freq_mean/SSize, 
             ErrorRange = sqrt(L1WidthAggregated$Freq_var/SSize^2 /
                                 L1WidthAggregated$Freq_N),
             TipWidth = 20)
lines(LengthVals, ExpL1WidthTa)
lines(LengthVals, ExpL1Width_nonTa)
par(new = T)
plot(LengthVals, SVals, type = "l", xaxt = "n", yaxt = "n", ylab = "", xlab = "",
     lty = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Selection coefficient')

CreateDisplayPdf('D:/L1polymORF/Figures/FreqVsL1Width_Group.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 4, width = 5)

# Plot individual points (length and frequency values)
par( mfrow = c(1, 1), oma = c(0.2,  0.2,  0.2,  0.2), 
     mai = c(1, 1, 0.2, 1),
     cex.lab = 1)
plot(L1TotData$L1width, 
     L1TotData$Freq/SSize, xlab = "LINE-1 length [bp]",
     ylab = "Mean LINE-1 frequency", col = rgb(0, 0, 0, alpha = 0.2), 
     pch = 16, ylim = c(0, 0.04))
t.test(Freq ~ blnFull, data = L1TotData)
L1FreqLengthSmoothed <- supsmu(L1TotData$L1width, 
                               L1TotData$Freq/SSize, bass = 5)
lines(L1FreqLengthSmoothed$x, L1FreqLengthSmoothed$y, lwd = 2, col = "green")
lines(LengthVals, ExpL1WidthTa, lwd = 2, col = "red")
lines(LengthVals, ExpL1Width_nonTa, lwd = 2, col = "red")
par(new = T)
plot(LengthVals, SVals, type = "l", xaxt = "n", yaxt = "n", ylab = "", 
     xlab = "", lwd = 2,col = "blue")
axis(side = 4)
mtext(side = 4, line = 3, 'Selection coefficient')

CreateDisplayPdf('D:/L1polymORF/Figures/FreqVsL1Width_smoothed.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 4, width = 5)

# Plot variance against width
VarL1Width <- sapply(SVals, function(x) VarAlleleFreq(x, N = 10^4, 
                                                      SampleSize = SSize,
                                                      DetectProb = 0.9,
                                                      LogRegCoeff = LogRegL1Ref$coefficients))
par( mfrow = c(1, 1), oma = c( 0.2,  0.2,  0.2,  0.2), 
     mai = c(1, 1, 0.2, 1),
     cex.lab = 1)
plot(L1WidthAggregated$L1width, 
     L1WidthAggregated$Freq_var / SSize^2, xlab = "LINE-1 length [bp]",
     ylab = "Variance LINE-1 frequency", ylim = c(0, 0.01))
lines(LengthVals, VarL1Width)
par(new = T)
plot(LengthVals, SVals, type = "l", xaxt = "n", yaxt = "n", ylab = "", xlab = "",
     lty = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Selection coefficient')

CreateDisplayPdf('D:/L1polymORF/Figures/FreqVarVsL1Width_Group.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 5, width = 5)


###############################################################
#                                                             #
#   Plot expected and observed allele frequency distribution  #
#                                                             #
###############################################################

FreqCountV <- 1:30

# Get per L1 width bin the site-frequency spectrum (SFS) and calculate a 
# multinomial probability per SFS
LBinsUnique <- unique(L1TotData$InsLengthClass)
LBinsUnique <- LBinsUnique[!is.na(LBinsUnique)]
SFS_list <- lapply(LBinsUnique, function(x){
  blnX   <- L1TotData$InsLengthClass == x
  Counts <- L1TotData$Freq[blnX]
  Counts <- Counts[!is.na(Counts)]
  if (length(Counts) > 0){
    SFS <- sapply(FreqCountV, function(y) sum(Counts == y))
  }
})
names(SFS_list) <- LBinsUnique

# Match data fram to SFS_list
binMatch <- match(L1WidthAggregated$InsLengthClass, names(SFS_list))
SFS_list <- SFS_list[binMatch]
length(SFS_list)
# Create a matrix of log probabilities of combinations of population 
# frequencies (rows) and selection coefficients (columns)  
cat("Calculating log probabilities of combinations of selection coefficients\n",
 "and population frequencies ...")
ExpMat <- sapply(1:nrow(L1WidthAggregated), function(x){
  L1WidthAggregated$Freq_N[x]*
  sapply(FreqCountV, function(y) {
    exp(AlleleFreqSample_pracma(k = y, s = L1WidthAggregated$sVals_mean[x], N = PopSize,
                     SampleSize = SSize,
                     DetectProb = 0.8,
                     LogRegCoeff = LogRegL1Ref$coefficients,
                     blnIns = T))
  })
})
cat("done!\n")


# Sum probabilities per frequency value to get expected frequencies
ExpFreq <- rowSums(ExpMat, na.rm = T)
plot(ExpFreq)


# Plot expected and observed frequency per length class
par(mfrow = c(4, 4), oma = c(5,  5,  0.2,  0.2), 
    mai = c(0.1, 0.1, 0.5, 0.2), cex.main = 0.75)

# Plot histogram of frequencies for all width class 
HF <- hist(L1_1000G$Frequency * SSize, breaks = 0:5000,
           xlim = c(0, 30), main = "All", xlab = "",
           ylab = "")
lines(HF$mids[1:30], ExpFreq)
for (i in 1:nrow(L1WidthAggregated)){
  x <- L1WidthAggregated$InsLengthClass[i]
  blnX   <- L1TotData$InsLengthClass == x
  Counts <- L1TotData$Freq[blnX]
  HF <- hist(Counts, breaks = 0:5000,
             xlim = c(0, 30), main = x, xlab = "",
             ylab = "")
  lines(HF$mids[1:30], ExpMat[,i])
  
}
mtext(side = 1, line = 3, 'Population frequency', outer = T)
mtext(side = 2, line = 3, 'Number of LINE-1s', outer = T)

CreateDisplayPdf('D:/L1ManuscriptFigures/ObsExpFreq.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 5, width = 5)

###############################################################
#                                                             #
#   Plot expected and observed allele frequency distribution  #
#   parameters esytimated insertion only  #
#                                                             #
###############################################################

FreqCountV <- 1:30

# Get per L1 width bin the site-frequency spectrum (SFS) and calculate a 
# multinomial probability per SFS
LBinsUnique <- unique(L1TotData$InsLengthClass)
LBinsUnique <- LBinsUnique[!is.na(LBinsUnique)]
SFS_list <- lapply(LBinsUnique, function(x){
  blnX   <- L1TotData$InsLengthClass == x
  Counts <- L1TotData$Freq[blnX]
  Counts <- Counts[!is.na(Counts)]
  if (length(Counts) > 0){
    SFS <- sapply(FreqCountV, function(y) sum(Counts == y))
  }
})
names(SFS_list) <- LBinsUnique

# Match data fram to SFS_list
binMatch <- match(L1WidthAggregated$InsLengthClass, names(SFS_list))
SFS_list <- SFS_list[binMatch]
length(SFS_list)
par(mfrow = c(1, 1))
plot(L1WidthAggregated$sVals2_mean, L1WidthAggregated$sVals_mean)
lines(c(0, 1), c(0, 1))
# Create a matrix of log probabilities of combinations of population 
# frequencies (rows) and selection coefficients (columns)  
cat("Calculating log probabilities of combinations of selection coefficients\n",
    "and population frequencies ...")
ExpMat <- sapply(1:nrow(L1WidthAggregated), function(x){
  L1WidthAggregated$Freq_N[x]*
    sapply(FreqCountV, function(y) {
      exp(AlleleFreqSample_pracma(k = y, s = L1WidthAggregated$sVals2_mean[x], N = PopSize,
                                  SampleSize = SSize,
                                  DetectProb = 0.8,
                                  LogRegCoeff = LogRegL1Ref$coefficients,
                                  blnIns = T))
    })
})
cat("done!\n")


# Sum probabilities per frequency value to get expected frequencies
ExpFreq <- rowSums(ExpMat, na.rm = T)
plot(ExpFreq)


# Plot expected and observed frequency per length class
par(mfrow = c(4, 4), oma = c(5,  5,  0.2,  0.2), 
    mai = c(0.1, 0.1, 0.5, 0.2), cex.main = 0.75)

# Plot histogram of frequencies for all width class 
HF <- hist(L1_1000G$Frequency * SSize, breaks = 0:5000,
           xlim = c(0, 30), main = "All", xlab = "",
           ylab = "")
lines(HF$mids[1:30], ExpFreq)
for (i in 1:nrow(L1WidthAggregated)){
  x <- L1WidthAggregated$InsLengthClass[i]
  blnX   <- L1TotData$InsLengthClass == x
  Counts <- L1TotData$Freq[blnX]
  HF <- hist(Counts, breaks = 0:5000,
             xlim = c(0, 30), main = x, xlab = "",
             ylab = "")
  lines(HF$mids[1:30], ExpMat[,i])
  
}
mtext(side = 1, line = 3, 'Population frequency', outer = T)
mtext(side = 2, line = 3, 'Number of LINE-1s', outer = T)

CreateDisplayPdf('D:/L1ManuscriptFigures/ObsExpFreq2.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 5, width = 5)

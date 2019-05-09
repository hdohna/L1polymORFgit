# The script below estimates selection coefficients of L1 from the 
# 1000 genome data

##########################################
#                                        #
#     Load packages                      #
#                                        #
##########################################

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(pracma)

##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

# Specify file paths
InputPath <- "D:/L1polymORF/Data/L1SelectionResults.RData"
InputPath <- "D:/L1polymORF/Data/L1SelectionResults_MELT.RData"
InputPath <- "D:/L1polymORF/Data/L1SelectionResults_MELT_Single.RData"

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("\n\nLoading and processing data ...")

# Load previously generated objects
load(InputPath)

###################################################
#                                                 #
#   Plot density vs. selection coefficient        #
#                                                 #
###################################################


# # Create a vector of selection coefficients
# SCoeffVect <- c(Promoter = ML_L1ExonIntron$par[1],
#                 Exon = sum(ML_L1ExonIntron$par[c(1, 2)]),
#                 Intron = sum(ML_L1ExonIntron$par[c(1, 3)]),
#                 Intergenic = ML_L1ExonIntron$par[1])
# names(SCoeffVect) <- sapply(names(SCoeffVect), 
#                             function(x) strsplit(x, "\\.")[[1]][1])
# 
# # Plot selection coefficient against 
# if (!all(names(SCoeffVect) == colnames(InsPerbp))){
#   stop("Selection coefficients and L1 densities are not in same order!")
# }
# if (!all(names(SCoeffVect) == names(MeanFreqs))){
#   stop("Selection coefficients and L1 frequencies are not in same order!")
# }
# 
# # Get sample size and create a range of s-values
# SSize <- 2*2504
# SVals <- seq(-0.0025, -0.00001, 0.00001)
# 
# # Probability of inclusion as funtion of selection coefficient
# PInclFun <- function(s, N = 10^4, N1, Nnf, SampleSize){
#   Pf1 <- ProbFix1(s, N = N) # Probability of fixation at 1
#   Pf  <- N1 / (Pf1*Nnf + N1) # Probability of fixation
#   # Probability of inclusion | no fixation
#   ProbL1 <- ProbAlleleIncluded(s, N = N, SampleSize = SampleSize) 
#   (1 - Pf)*ProbL1 #+ Pf * Pf1
#   
# }
# PIncl <- sapply(SVals, function(s) PInclFun(s, N1 = N1, Nnf = Nnf, SampleSize = SSize)) # Probability of fixation at 1
# ExpL1 <- sapply(SVals, function(x) ExpAlleleFreq(x, N = 10^4, SampleSize = SSize))
# 
# # Get probability of inclusion and mean frequency as function of selection
# # coefficient for each insertion region
# PInclVect <- sapply(SCoeffVect, function(s) {
#   PInclFun(s, N1 = N1, Nnf = Nnf, SampleSize = SSize)})
# SqDiffDens <- function(x) {sum((InsPerbp[2,] - x *PInclVect)^2)}
# CoeffVals <- seq(10^5, 5*10^5, 100)
# plot(CoeffVals, sapply(CoeffVals, function(x) SqDiffDens(x)), type = "l")
# OptCoeff <- optim(par = 3 * mean(InsPerbp[2,]) / mean(PIncl), 
#                   fn = function(x) SqDiffDens(x),
#                   method = "Brent",
#                   lower = 10^4, upper = 10^6)
# 
# plot(InsPerbp[2,], PInclVect* 3 * mean(InsPerbp[2,]) / mean(PIncl))
# lines(c(0, 10), c(0, 10))
# 
# par(oma = c(2, 1, 1, 3), mfcol = c(2, 2), mai = c(1, 1, 0.2, 0.2),
#     cex.axis = 1, cex.lab = 1.5)
# layout(rbind(1:2, c(3, 3)), widths = c(1, 1))
# layout(rbind(c(1, 1, 2, 2), c(0, 3, 3, 0)), widths = c(1, 1))
# 
# # Plot LINE-1 frequency against number of LINE-1 per Mb
# plot(MeanFreqs, InsPerbp[2,], ylab = "LINE-1s per Mb", 
#      xlab = "Mean LINE-1 frequency", main = "A", ylim = c(0, 3))
# lines(ExpL1, PIncl * OptCoeff$par)
# text(MeanFreqs + c(3*10^-3, 3*10^-3, 3*10^-3, -3*10^-3),
#      InsPerbp[2,] + 2*10^(-1)*c(1, 0, 0, -1.2), 
#      names(SCoeffVect))
# 
# # Plot expected frequency versus observed mean frequency
# plot(SCoeffVect, MeanFreqs, ylab = "Mean LINE-1 frequency", 
#      xlab = "Selection coefficient", xlim = c(-0.0025, 0.0007), main = "B")
# text(SCoeffVect + 2*c(0.0003, 0, -0.0003, -0.0003), 
#      MeanFreqs + c(0, 3*10^-3, 0, 0), names(SCoeffVect))
# lines(SVals, ExpL1)
# 
# 
# # Plot probability for inclusion versus number of LINE-1 per Mb
# plot(SCoeffVect, InsPerbp[2,], ylab = "LINE-1s per Mb", 
#      xlab = "Selection coefficient", xlim = c(-0.0025, 0), ylim = c(0, 3),
#      main = "C")
# text(SCoeffVect, InsPerbp[2,] + 2*10^(-1), names(SCoeffVect))
# par(new = TRUE)
# plot(SVals, PIncl, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",
#      ylim = c(0, 5*10^-6))
# axis(side = 4)
# mtext("Inclusion probability", 4, line = 3)
# #mtext("Selection coefficient", 1, line = 3)
# 
# CreateDisplayPdf('D:/L1polymORF/Figures/SelectionPerRegion_MELT.pdf',
#                  PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
#                  height = 7, width = 7)

###################################################
#                                                 #
#   Plot frequency vs. insertion length           #
#                                                 #
###################################################

# Create a vector of L1 start classes
L1TotData$InsLengthClass <- cut(L1TotData$L1width, breaks = 
                                  seq(0, 6750, 750))
min(L1TotData$Freq, na.rm = T)
blnNoZero <- L1TotData$Freq > 0
# Get mean L1 frequency per start
L1WidthAggregated <- aggregate(L1TotData[blnNoZero,c("L1width", "Freq")], 
                               by = list(L1TotData$InsLengthClass[blnNoZero]), 
                               FUN = function(x) mean(x, na.rm = T))
L1WidthAggregated_var <- aggregate(L1TotData[blnNoZero,c("L1width", "Freq")], 
                                   by = list(L1TotData$InsLengthClass[blnNoZero]), 
                                   FUN = function(x) var(x, na.rm = T))
L1WidthAggregated_n <- aggregate(L1TotData[blnNoZero,c("L1width", "Freq")], 
                                 by = list(L1TotData$InsLengthClass[blnNoZero]), 
                                 FUN = function(x) sum(!is.na(x)))
L1WidthAggregated <- aggregate(L1TotData[L1TotData$blnIns & blnNoZero,c("L1width", "Freq")],
                               by = list(L1TotData$InsLengthClass[L1TotData$blnIns & blnNoZero]),
                               FUN = function(x) mean(x, na.rm = T))
L1WidthAggregated_var <- aggregate(L1TotData[L1TotData$blnIns & blnNoZero,c("L1width", "Freq")],
                               by = list(L1TotData$InsLengthClass[L1TotData$blnIns & blnNoZero]),
                               FUN = function(x) var(x, na.rm = T))
L1WidthAggregated_n <- aggregate(L1TotData[L1TotData$blnIns & blnNoZero,c("L1width", "Freq")],
                                   by = list(L1TotData$InsLengthClass[L1TotData$blnIns & blnNoZero]),
                                   FUN = function(x) sum(!is.na(x)))
# 
# Get sample size and create a range of s-values
SSize <- 2*2504
StartVals  <- seq(0, 6200, 200)
Full       <- StartVals >= 6000
SVals <- ML_L1widthL1full$par[1] + ML_L1widthL1full$par[2]*StartVals +
  ML_L1widthL1full$par[3]*Full
DetectProb <- exp(L1SizeDetectCoeff[1] + 
      L1SizeDetectCoeff[2] * StartVals) / 
  (1 + exp(L1SizeDetectCoeff[1] + 
             L1SizeDetectCoeff[2] * StartVals))

SVals2 <- -0.000001 - 2*10^-7*StartVals +
  1.08*10^-3*Full

# Calculate expected frequency per 
ExpL1Width <- sapply(1:length(SVals), function(i) ExpAlleleFreq(s = SVals[i], N = PopSize, 
                                                      SampleSize = MEInsSamplesize,
                                                      DetectProb = DetectProb[i],
                                                      blnIns = T, 
                                                      LogRegCoeff = LogRegL1Ref$coefficients))
# ExpL1Width2 <- sapply(SVals2, function(x) ExpAlleleFreq(s = x, N = 10^5, 
#                                                       SampleSize = 2*2504,
#                                                       DetectProb = 0.9,
#                                                       LogRegCoeff = LogRegL1Ref$coefficients))
par( mfrow = c(1, 1), oma = c( 0.2,  0.2,  0.2,  0.2), 
     mai = c(1, 1, 0.2, 1),
     cex.lab = 1)
plot(L1WidthAggregated$L1width, 
     L1WidthAggregated$Freq/SSize, xlab = "LINE-1 length [bp]",
     ylab = "Mean LINE-1 frequency", ylim = c(0.0005, 0.004))
AddErrorBars(MidX = L1WidthAggregated$L1width, 
             MidY = L1WidthAggregated$Freq/SSize, 
             ErrorRange = sqrt(L1WidthAggregated_var$Freq/SSize^2 /
                                 L1WidthAggregated_n$Freq),
             TipWidth = 20)
#lines(StartVals, ExpL1Width2)
lines(StartVals, ExpL1Width)
par(new = T)
plot(StartVals, SVals, type = "l", xaxt = "n", yaxt = "n", ylab = "", xlab = "",
     lty = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Selection coefficient')

CreateDisplayPdf('D:/L1polymORF/Figures/FreqVsL1Width.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 4, width = 5)

# Plot Individual points
par( mfrow = c(1, 1), oma = c( 0.2,  0.2,  0.2,  0.2), 
     mai = c(1, 1, 0.2, 1),
     cex.lab = 1)
plot(L1TotData$L1width, 
     L1TotData$Freq/SSize, xlab = "LINE-1 length [bp]",
     ylab = "Mean LINE-1 frequency", col = rgb(0, 0, 0, alpha = 0.2), 
     pch = 16)
L1FreqLengthSmoothed <- supsmu(L1TotData$L1width, 
                               L1TotData$Freq/SSize)
lines(L1FreqLengthSmoothed$x, L1FreqLengthSmoothed$y, col = "red")
lines(StartVals, ExpL1Width, lty = 4, lwd = 2, col = "red")
par(new = T)
plot(StartVals, SVals, type = "l", xaxt = "n", yaxt = "n", ylab = "", xlab = "",
     col = "blue")
axis(side = 4, col = "blue")
mtext(side = 4, line = 4, 'Selection coefficient')

CreateDisplayPdf('D:/L1polymORF/Figures/FreqVsL1Width_smoothed.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 6, width = 6)

# Plot variance against width
VarL1Width <- sapply(SVals, function(x) VarAlleleFreq(x, N = 10^4, 
                                                      SampleSize = 2*2504,
                                                      DetectProb = 0.9,
                                                      LogRegCoeff = LogRegL1Ref$coefficients))
par( mfrow = c(1, 1), oma = c( 0.2,  0.2,  0.2,  0.2), 
     mai = c(1, 1, 0.2, 1),
     cex.lab = 1)
plot(L1WidthAggregated$InsLength, 
     L1WidthAggregated_var$Frequency, xlab = "LINE-1 length [bp]",
     ylab = "Variance LINE-1 frequency")
lines(StartVals, VarL1Width)
par(new = T)
plot(StartVals, SVals, type = "l", xaxt = "n", yaxt = "n", ylab = "", xlab = "",
     lty = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Selection coefficient')

CreateDisplayPdf('D:/L1polymORF/Figures/FreqVsL1Width.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 5, width = 5)


###############################################################
#                                                             #
#   Plot expected and observed allele frequency distribution  #
#                                                             #
###############################################################

# Create a matrix of predictor variables (L1 width and boolean variable for
# for full-length L1)
XMat <- as.matrix(cbind(1, L1TotData[, c("L1width", "blnFull")]))
blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(XMat[x,])))


sVals <- XMat %*% ML_L1widthL1full$par
sVals <- sVals[!is.na(sVals)]

freqVals <- seq(1, SSize, 100)

ExpFreqMat <- sapply(sVals, function(x){
  sapply(freqVals, function(y) {
    AlleleFreqSample(k = y, s = x, N = 10^4,
                     SampleSize = SSize, 
                     DetectProb = 0.9,
                     LogRegCoeff = LogRegL1Ref$coefficients, 
                     blnIns = T)
  })
})
dim(ExpFreqMat)
ExpFreq <- rowSums(exp(ExpFreqMat))
LogFreq <- log(ExpFreq)
plot(ExpFreq, type = "l")
plot(LogFreq, type = "l")
plot(LogFreq, type = "l", ylim = c(-10, 5))

HF <- hist(L1_1000G$Frequency)
plot(HF$counts)
plot(log(HF$counts))

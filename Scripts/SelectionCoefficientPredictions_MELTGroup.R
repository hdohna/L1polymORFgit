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
#library(pracma)

##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

# Specify file paths
InputPath <- "D:/L1polymORF/Data/L1SelectionResults_MELT_GroupwithSim.RData"


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
#   Plot frequency vs. insertion length           #
#                                                 #
###################################################

# Create a vector of L1 start classes
L1TotData$InsLengthClass <- cut(L1TotData$L1width, breaks = 
                                  seq(0, 6500, 5000))

# Get mean L1 frequency per start
L1WidthAggregated <- AggDataFrame(L1TotData, 
                                  GroupCol = "InsLengthClass", 
                                  MeanCols = c("L1width", "Freq"), 
                                  LengthCols = "Freq",
                                  VarCols = "Freq")
                                              
# Get sample size and create a range of s-values
SSize      <- 2*2504
StartVals  <- seq(0, 6200, 200)
Full       <- StartVals >= 6000
ModelFit1
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
     ylab = "Mean LINE-1 frequency", ylim = c(0.0005, 0.009))
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
t.test(Freq ~ blnFull, data = L1TotData)
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
plot(L1WidthAggregated$L1width, 
     L1WidthAggregated_var$Freq, xlab = "LINE-1 length [bp]",
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

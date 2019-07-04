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

# TODO: Find out why expected frequency value based on estimated selection 
# coefficient is always higher than observed frequency value

# Create a matrix of predictor variables (L1 width and boolean variable for
# for full-length L1)
XMat  <- as.matrix(cbind(1, L1TotData[, c("L1width", "blnFull")]))

# Create a vector of selection coefficients per L1
L1TotData$sVals <- XMat %*% ModelFit1$ML_abc$par

# Get sample size and create a range of s-values
SSize      <- L1TotData$SampleSize[1]
LengthVals <- seq(0, 6200, 200)
Full      <- LengthVals >= 6000
SVals <- ModelFit1$ML_abc$par[1] + ModelFit1$ML_abc$pa[2]*Full +
  ModelFit1$ML_abc$par[3]*LengthVals
DetectProb <- L1TotData$DetectProb[1]

# Calculate expected frequency per L1 width
ExpL1Width <- sapply(1:length(SVals), function(i) {
  ExpAlleleFreq_pracma(s = SVals[i], N = PopSize, SampleSize = MEInsSamplesize,
                DetectProb = DetectProb, blnIns = T, 
                LogRegCoeff = LogRegL1Ref$coefficients)
  })

# Create a vector of L1 length classes
L1TotData$InsLengthClass <- cut(L1TotData$L1width, breaks = 
                                  seq(0, 6750, 750))

# Get mean L1 frequency per length
L1WidthAggregated <- AggDataFrame(L1TotData, 
                                  GroupCol = "InsLengthClass", 
                                  MeanCols = c("L1width", "Freq", "blnFull", "sVals"), 
                                  LengthCols = "Freq",
                                  VarCols = "Freq")
idxL1LengthMatch <- match(L1TotData$InsLengthClass, 
                          L1WidthAggregated$InsLengthClass)                                              
  
# Plot mean frequency, expected frequency
par(mfrow = c(1, 1), oma = c( 0.2,  0.2,  0.2,  0.2), 
     mai = c(1, 1, 0.2, 1.5),
     cex.lab = 1)
plot(L1WidthAggregated$L1width_mean, 
     L1WidthAggregated$Freq_mean/SSize, xlab = "LINE-1 length [bp]",
     ylab = "Mean LINE-1 frequency", ylim = c(0, 0.03))
AddErrorBars(MidX = L1WidthAggregated$L1width_mean, 
             MidY = L1WidthAggregated$Freq_mean/SSize, 
             ErrorRange = sqrt(L1WidthAggregated$Freq_var/SSize^2 /
                                 L1WidthAggregated$Freq_N),
             TipWidth = 20)
lines(LengthVals, ExpL1Width)
par(new = T)
plot(LengthVals, SVals, type = "l", xaxt = "n", yaxt = "n", ylab = "", xlab = "",
     lty = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Selection coefficient')

CreateDisplayPdf('D:/L1polymORF/Figures/FreqVsL1Width_Group.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 4, width = 5)

# Plot individual points (length and frequency values)
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
lines(LengthVals, ExpL1Width, lty = 4, lwd = 2, col = "red")
par(new = T)
plot(LengthVals, SVals, type = "l", xaxt = "n", yaxt = "n", ylab = "", xlab = "",
     col = "blue")
axis(side = 4, col = "blue")
mtext(side = 4, line = 4, 'Selection coefficient')

CreateDisplayPdf('D:/L1polymORF/Figures/FreqVsL1Width_smoothed.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 6, width = 6)

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


# Create a vector of population frequencies
freqVals <- seq(1, SSize + 1, 100)

# Create a matrix of log probabilities of combinations of population 
# frequencies (rows) and selection coefficients (columns)  
cat("Calculating log probabilities of combinations of selection coefficients\n",
 "and population frequencies ...")
LogProbMat <- sapply(1:nrow(L1WidthAggregated), function(x){
  L1WidthAggregated$Freq_N[x]*
  sapply(1:50, function(y) {
    exp(AlleleFreqSample(k = y, s = L1WidthAggregated$sVals_mean[x], N = PopSize,
                     SampleSize = SSize,
                     DetectProb = 0.8,
                     LogRegCoeff = LogRegL1Ref$coefficients,
                     blnIns = T))
  })
})
cat("done!\n")


# Sum probabilities per frequency value to get expected frequencies
ExpFreq <- rowSums(LogProbMat, na.rm = T)
plot(ExpFreq)

# Plot histogram of frequencies 
freqValsHist <- freqVals
freqValsHist[1] <- 0
HF <- hist(L1_1000G$Frequency * SSize, breaks = 0:5000,
           xlim = c(0, 50))
lines(freqVals[-length(freqVals)], ExpFreq / sum(ExpFreq)/SSize * 50)

plot(HF$counts)
plot(log(HF$counts))

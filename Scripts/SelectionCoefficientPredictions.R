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


# Create a vector of selection coefficients
SCoeffVect <- c(Promoter = ML_L1ExonIntron$par[1],
                Exon = sum(ML_L1ExonIntron$par[c(1, 2)]),
                Intron = sum(ML_L1ExonIntron$par[c(1, 3)]),
                Intergenic = ML_L1ExonIntron$par[1])
names(SCoeffVect) <- sapply(names(SCoeffVect), 
                            function(x) strsplit(x, "\\.")[[1]][1])

# Plot selection coefficient against 
if (!all(names(SCoeffVect) == colnames(InsPerbp))){
  stop("Selection coefficients and L1 densities are not in same order!")
}
if (!all(names(SCoeffVect) == names(MeanFreqs))){
  stop("Selection coefficients and L1 frequencies are not in same order!")
}

# Get sample size and create a range of s-values
SSize <- 2*2504
SVals <- seq(-0.0025, -0.00001, 0.00001)

# Probability of inclusion as funtion of selection coefficient
PInclFun <- function(s, N = 10^4, N1, Nnf, SampleSize){
  Pf1 <- ProbFix1(s, N = N) # Probability of fixation at 1
  Pf  <- N1 / (Pf1*Nnf + N1) # Probability of fixation
  # Probability of inclusion | no fixation
  ProbL1 <- ProbAlleleIncluded(s, N = N, SampleSize = SampleSize) 
  (1 - Pf)*ProbL1 #+ Pf * Pf1
  
}
PIncl <- sapply(SVals, function(s) PInclFun(s, N1 = N1, Nnf = Nnf, SampleSize = SSize)) # Probability of fixation at 1
ExpL1 <- sapply(SVals, function(x) ExpAlleleFreq(x, N = 10^4, SampleSize = SSize))

# Get probability of inclusion and mean frequency as function of selection
# coefficient for each insertion region
PInclVect <- sapply(SCoeffVect, function(s) {
  PInclFun(s, N1 = N1, Nnf = Nnf, SampleSize = SSize)})
SqDiffDens <- function(x) {sum((InsPerbp[2,] - x *PInclVect)^2)}
CoeffVals <- seq(10^5, 5*10^5, 100)
plot(CoeffVals, sapply(CoeffVals, function(x) SqDiffDens(x)), type = "l")
OptCoeff <- optim(par = 3 * mean(InsPerbp[2,]) / mean(PIncl), 
                  fn = function(x) SqDiffDens(x),
                  method = "Brent",
                  lower = 10^4, upper = 10^6)

plot(InsPerbp[2,], PInclVect* 3 * mean(InsPerbp[2,]) / mean(PIncl))
lines(c(0, 10), c(0, 10))

par(oma = c(2, 1, 1, 3), mfcol = c(2, 2), mai = c(1, 1, 0.2, 0.2),
    cex.axis = 1, cex.lab = 1.5)
layout(rbind(1:2, c(3, 3)), widths = c(1, 1))
layout(rbind(c(1, 1, 2, 2), c(0, 3, 3, 0)), widths = c(1, 1))

# Plot LINE-1 frequency against number of LINE-1 per Mb
plot(MeanFreqs, InsPerbp[2,], ylab = "LINE-1s per Mb", 
     xlab = "Mean LINE-1 frequency", main = "A", ylim = c(0, 3))
lines(ExpL1, PIncl * OptCoeff$par)
text(MeanFreqs + c(3*10^-3, 3*10^-3, 3*10^-3, -3*10^-3),
     InsPerbp[2,] + 2*10^(-1)*c(1, 0, 0, -1.2), 
     names(SCoeffVect))

# Plot expected frequency versus observed mean frequency
plot(SCoeffVect, MeanFreqs, ylab = "Mean LINE-1 frequency", 
     xlab = "Selection coefficient", xlim = c(-0.0025, 0.0007), main = "B")
text(SCoeffVect + 2*c(0.0003, 0, -0.0003, -0.0003), 
     MeanFreqs + c(0, 3*10^-3, 0, 0), names(SCoeffVect))
lines(SVals, ExpL1)


# Plot probability for inclusion versus number of LINE-1 per Mb
plot(SCoeffVect, InsPerbp[2,], ylab = "LINE-1s per Mb", 
     xlab = "Selection coefficient", xlim = c(-0.0025, 0), ylim = c(0, 3),
     main = "C")
text(SCoeffVect, InsPerbp[2,] + 2*10^(-1), names(SCoeffVect))
par(new = TRUE)
plot(SVals, PIncl, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",
     ylim = c(0, 5*10^-6))
axis(side = 4)
mtext("Inclusion probability", 4, line = 3)
#mtext("Selection coefficient", 1, line = 3)

CreateDisplayPdf('D:/L1polymORF/Figures/SelectionPerRegion_MELT.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7, width = 7)

###################################################
#                                                 #
#   Plot frequency vs. insertion length           #
#                                                 #
###################################################

# Create a vector of L1 start classes
L1_1000G$InsLengthClass <- cut(L1_1000G$InsLength, breaks = 
                                  seq(0, 6500, 500))

# Get mean L1 frequency per start
L1WidthAggregated <- aggregate(L1_1000G[,c("InsLength", "Frequency")], 
                               by = list(L1_1000G$InsLengthClass), FUN = mean)
L1WidthAggregated_var <- aggregate(L1_1000G[,c("InsLength", "Frequency")], 
                               by = list(L1_1000G$InsLengthClass), FUN = var)
L1WidthAggregated_n <- aggregate(L1_1000G[,c("InsLength", "Frequency")], 
                                   by = list(L1_1000G$InsLengthClass), FUN = length)

# Get sample size and create a range of s-values
SSize <- 2*2504
StartVals  <- seq(0, 6100, 10)
Full       <- StartVals >= 6000
SVals <- ML_L1widthL1full$par[1] + ML_L1widthL1full$par[2]*StartVals +
  ML_L1widthL1full$par[3]*Full

# Plot expected frequency versus observed mean frequency
ExpL1Width <- sapply(SVals, function(x) ExpAlleleFreq(x, N = 10^4, 
                                                      SampleSize = 2*2504))
par( mfrow = c(1, 1))
plot(L1WidthAggregated$InsLength, 
     L1WidthAggregated$Frequency, xlab = "LINE-1 length [bp]",
     ylab = "Mean LINE-1 frequency", ylim = c(0, 0.05))
AddErrorBars(MidX = L1WidthAggregated$InsLength, 
             MidY = L1WidthAggregated$Frequency, 
             ErrorRange = sqrt(L1WidthAggregated_var$Frequency /
                                 L1WidthAggregated_n$Frequency),
             TipWidth = 20)
lines(StartVals, ExpL1Width)
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


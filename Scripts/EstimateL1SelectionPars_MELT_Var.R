# The script below estimates selection coefficients of L1 from the 
# 1000 genome data using insertion estimates obtained by MELT 
# 
# Load previous results
load("D:/L1polymORF/Data/L1SelectionResults_MELT.RData")
library(pracma)

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Estimate maximum likelihood for a single selection coefficient
ML_1ParVar <-  constrOptim(theta = c(a = ML_1Par$par, SD = 10^-3),
                           f = function(x) -AlleleFreqLogLik_4Par(
                             Freqs = round(L1TotData$Freq[!blnNA], 0),
                             Counts = rep(1, sum(!blnNA)),
                             Predict = PredictMat[!blnNA, 1:3],
                             a = x[1], b = 0, c = 0, d = 0, SD = x[2], N = 10^4,
                             SampleSize = L1TotData$SampleSize[!blnNA],
                             blnIns = L1TotData$blnIns[!blnNA], 
                             DetectProb = 0.9, VariableSelection = T,
                             LowerS = -0.5,
                             UpperS = 0.5),
                           grad = NULL,
                           ui = rbind(c(1, 0),  c(0, 1),   
                                      c(-1, 0), c(0, -1)),
                           ci = c(a = -0.02, SD = 5*10^-4, 
                                  a = -0.02, SD = -1),
                           method = "Nelder-Mead")
cat("done!\n")

Freqs = round(L1TotData$Freq[!blnNA], 0)
Counts = rep(1, sum(!blnNA))
Predict = PredictMat[!blnNA, 1:3]
a = ML_1Par$par
b = 0
c = 0
d = 0
SD = 5*10^-4
N = 10^4
SampleSize = L1TotData$SampleSize[!blnNA]
blnIns = L1TotData$blnIns[!blnNA]
DetectProb = 0.9 
VariableSelection = T
i <- 1

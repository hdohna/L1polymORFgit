# The script below estimates selection coefficients of L1 from the 
# 1000 genome data using insertion estimates obtained by MELT 
# 
# Load previous results
load("D:/L1polymORF/Data/L1SelectionResults_MELT.RData")
library(pracma)

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Estimate maximum likelihood for a single selection coefficient and its 
# variance
cat("Maximizing likelihood for one selection parameter and variance ...")
MML_1ParVar <-  constrOptim(theta = c(a = ML_1Par$par, SD = 10^-3),
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
# Best estimate so far a = -0.0004546913 SD = 0.001621152, log-lik = -13783.64

# Estimate maximum likelihood for three selection coefficient parameters 
# and selection coefficient variance
cat("Maximizing likelihood for three selection parameters and variance ...")
ML_L1startL1full_Var <- constrOptim(theta = c(a = MML_1ParVar$par[1], 
                                              c = ML_L1startL1full$par[2], 
                                              d = ML_L1startL1full$par[3], 
                                              SD = 0.9*10^-3),
                                f = function(x) -AlleleFreqLogLik_4Par(
                                  Freqs = round(L1TotData$Freq[!blnNA], 0),
                                  Counts = rep(1, sum(!blnNA)),
                                  Predict = PredictMat[!blnNA, 1:3],
                                  a = x[1], b = 0, c = x[2], d = x[3], SD = x[4], 
                                  N = 10^4, 
                                  SampleSize = L1TotData$SampleSize[!blnNA],
                                  blnIns = L1TotData$blnIns[!blnNA], 
                                  DetectProb = 0.9, VariableSelection = T,
                                  LowerS = -0.5,
                                  UpperS = 0.5),
                                grad = NULL,
                                ui = rbind(c(1, 0, 0, 0),  c(0, 1, 0, 0),  c(0, 0, 1, 0), c(0, 0, 0, 1),
                                           c(-1, 0, 0, 0), c(0, -1, 0, 0), c(0, 0, -1, 0), c(0, 0, 0, -1)),
                                ci = c(a = -0.01, b = -10^(-6), d = -10^(-3), SD = 5*10^-4,
                                       a = -0.02, b = -10^(-6), d = -10^(-3), SD = -1),
                                method = "Nelder-Mead")


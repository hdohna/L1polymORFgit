# The script below estimates selection coefficients of L1 from the 
# 1000 genome data USING SAMPLED L1 LENGTHS!

##########################################
#                                        #
#     Load packages                      #
#                                        #
##########################################

# Source start script
source('/home/hzudohna/L1polymORFgit/_Start_L1polymORF_scg4.R')

# Load packages
library(GenomicRanges)
library(pracma)

# Load data
load("/labs/dflev/hzudohna/1000Genomes/L1SelectionResults_MELT_GroupwithSim.RData")


# Fit model for various samples
ModelFitSamples_pracma <- lapply(1:100, function(x) {
  
  cat("\n********   Estimating effect of insertion length: sample", x, "   **********\n")
  blnL1widthNA <- is.na(L1TotData$L1width)
  L1TotData$L1width_sample <- NA
  L1TotData$L1width_sample[!blnL1widthNA] <- SampleTrueL1Width(
    SimL1widthTrue = L1DetectAgg_withL1$L1widthTrue, 
    SimL1widthEst  = L1DetectAgg_withL1$L1widthEst,
    SimL1Freq = L1DetectAgg_withL1$EstFreq,
    EstL1width = L1TotData$L1width[!blnL1widthNA],
    ObsL1Freq = L1TotData$Freq[!blnL1widthNA],
    L1widthBreaks = seq(0, 6500, 500))
   
  # Create a matrix of predictor variables (L1 start and boolean variable for)
  PredictMat <- L1TotData[, c("blnFull", "L1width_sample", "Freq")]
  blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))
  
  FitSelectionModels_pracma(PredictMat[!blnNA, ],  
                                  Freqs = round(L1TotData$Freq[!blnNA], 0), 
                                  Counts = rep(1, sum(!blnNA)), 
                                  PopSize = PopSize, 
                                  SampleSize = L1TotData$SampleSize[!blnNA],
                                  blnIns = L1TotData$blnIns[!blnNA], 
                                  LogRegCoeff = LogRegL1Ref$coefficients,
                                  DetectProb = L1TotData$DetectProb[!blnNA],
                                  aBorder = 0.003, 
                                  bBorder = 10^(-2), 
                                  cBorder = 10^(-5))
  
})


# Save everything
save.image("/labs/dflev/hzudohna/1000Genomes/L1SelectionResults_MELT_GroupwithSim100Saples.RData")

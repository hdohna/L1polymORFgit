library(fields)
#load("D:/L1polymORF/Data/SelectionParameterFit_1000G_maxF0.6_Acc_Quant2.RData")
#load("D:/L1polymORF/Data/SelectionParameterFit_1000G_maxF0.3_FineGrid.RData")
load("D:/L1polymORF/Data/SelectionParameterFit_Bossinot_KS_FineGrid.RData")

##############################
#                            #
#     Define functions       #
#                            #
##############################

# Plot statistic difference against parameter values
PlotResults <- function(ResultList){
  proPs <- unique(ResultList$proPsRep)
  par(mfrow = c(2, 2))
  # Plot difference vs alpha
  Cols <- rainbow(length(proPs))
  plot(aVals, DiffKS, xlab = "alpha")
  for (i in 1:length(Cols)){
    blnProps <- ResultList$proPsRep == proPs[i]
    points(ResultList$aVals[blnProps], ResultList$DiffMean[blnProps], col = Cols[i])
  }
  legend("bottomright", legend = proPs, col = Cols, pch = 1, cex = 0.5)
  
  # Plot difference vs selection coefficient
  aValsBasic <- unique(ResultList$aVals)
  Cols <- rainbow(length(aValsBasic))
  plot(ResultList$proPsRep, ResultList$DiffMean, xlab = "Mean fitness")
  for (i in 1:length(Cols)){
    blnA <- ResultList$aVals == aValsBasic[i]
    points(ResultList$proPsRep[blnA], ResultList$DiffMean[blnA], col = Cols[i])
  }
  legend("bottomright", legend = aValsBasic, col = Cols, pch = 1, cex = 0.5)
  
  # Plot minimum difference per alpha
#  plot(ResultList$aVals[idxMinDiffMean], ResultList$DiffMean[idxMinDiffMean], xlab = "alpha")
  
}

# Function to plot heatmap
PlotHeatmap <- function(ResultList){
  DiffMat <- matrix(ResultList$GridSummary, 
                      nrow = length(unique(ResultList$FitMeansRep)))
  image.plot(x = unique(ResultList$FitMeansRep),
             y = unique(ResultList$FitVarsRep),
             z = DiffMat)
  
}

# Function to get the best-fitting parameter values
ParsBestFit <- function(ResultList){
  idxMinDiff <- which.min(ResultList$DiffMean)
  c(BestMean = ResultList$FitMeansRep[idxMinDiff],
    BestVar  = ResultList$FitVarsRep[idxMinDiff],
    MinDiff  = ResultList$DiffMean[idxMinDiff])
}


##############################
#                            #
#     Plot results           #
#                            #
##############################

par(mfrow = c(2, 2))
PlotHeatmap(ResultListFull_1000G_hist)
PlotHeatmap(ResultListFragm_1000G_hist)

##############################
#                            #
#     Explore parameter          #
#                            #
##############################

# Get best-fitting mean and variance
BestFitPars_Full <- ParsBestFit(ResultListFull_1000G_hist)
BestFitPars_Fragm <- ParsBestFit(ResultListFragm_1000G_hist)

# Plot parameter vs deviation
blnBestMean <- ResultListFull_1000G_hist$FitMeansRep == 

# Plot the difference between observed and simulated frequency
plot(ResultListFull_1000G_hist$GridSummary[1, ] - ResultListFull_1000G_hist$ObsSummary[1])
plot(ResultListFull_1000G_hist$GridSummary[2, ] - ResultListFull_1000G_hist$ObsSummary[2])
plot(ResultListFull_1000G_hist$GridSummary[3, ] - ResultListFull_1000G_hist$ObsSummary[3])
hist(ResultListFull_1000G_hist$GridSummary)
# Generate an example frequency distribution for a reasonable set of parameters
# a1    <- 100
# fitn1 <- 0.86
# a2    <- 20
# fitn2 <- 0.86
# SampleFreq11 <- GenerateAlleleFreq(Gshape = a1, GSscale = 1/a1 * fitn1,
#                                   NrGen = 10^5)
# SampleFreq12 <- GenerateAlleleFreq(Gshape = a1, GSscale = 1/a1 * fitn1,
#                                   NrGen = 10^3)
# SampleFreq21 <- GenerateAlleleFreq(Gshape = a2, GSscale = 1/a2 * fitn2,
#                                    NrGen = 10^5)
# SampleFreq22 <- GenerateAlleleFreq(Gshape = a2, GSscale = 1/a2 * fitn2,
#                                    NrGen = 10^3)
# par(mfrow = c(2, 2))
# hist(SampleFreq11, breaks = seq(0, 1, 0.02))
# hist(SampleFreq12, breaks = seq(0, 1, 0.02))
# hist(SampleFreq21, breaks = seq(0, 1, 0.02))
# hist(SampleFreq22, breaks = seq(0, 1, 0.02))
# hist(FreqFull_1000G, breaks = seq(0, 1, 0.02))
# hist(FreqFragm_1000G, breaks = seq(0, 1, 0.02))
# 
# ExploreSelectionParameterGrid

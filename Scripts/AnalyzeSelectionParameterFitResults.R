library(fields)
load("D:/L1polymORF/Data/SelectionParameterFit_1000G_maxF0.6_Acc_Quant2.RData")

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
  DiffKSMat <- matrix(ResultList$DiffKS, 
                      nrow = length(unique(ResultList$proPsRep)))
  image.plot(x = unique(ResultList$proPsRep),
        y = unique(ResultList$aVals),
        z = DiffKSMat)

}
PlotHeatmap <- function(ResultList){
  DiffKSMat <- matrix(ResultList$DiffQuantMean, 
                      nrow = length(unique(ResultList$proPsRep)))
  image.plot(x = unique(ResultList$proPsRep),
             y = unique(ResultList$aVals),
             z = DiffKSMat)
  
}

par(mfrow = c(2, 2))
PlotHeatmap(ResultList2Full_1000G)
PlotHeatmap(ResultList2Fragm_1000G)

ResultList2Full_1000G$DiffQuantMean
ResultList2Full_1000G$DiffQuant[, 1:10]

# Generate an example frequency distribution for a reasonable set of parameters
a    <- 100
fitn <- 0.86
SampleFreq1 <- GenerateAlleleFreq(Gshape = a, GSscale = 1/a * fitn,
                                  NrGen = 10^5)
par(mfrow = c(2, 2))
hist(SampleFreq1, breaks = seq(0, 1, 0.02))
hist(FreqFull_1000G, breaks = seq(0, 1, 0.02))
hist(FreqFragm_1000G, breaks = seq(0, 1, 0.02))

library(fields)
load("D:/L1polymORF/Data/SelectionParameterFit_1000G_MaxFreq0.5.RData")

PlotResults <- function(ResultList){
  proPs <- unique(ResultList$proPsRep)
  par(mfrow = c(2, 2))
  # Plot difference vs alpha
  Cols <- rainbow(length(proPs))
  plot(aVals, DiffKS, xlab = "alpha")
  for (i in 1:length(Cols)){
    blnProps <- ResultList$proPsRep == proPs[i]
    points(ResultList$aVals[blnProps], ResultList$DiffKS[blnProps], col = Cols[i])
  }
  legend("bottomright", legend = proPs, col = Cols, pch = 1, cex = 0.5)
  
  # Plot difference vs selection coefficient
  aValsBasic <- unique(ResultList$aVals)
  Cols <- rainbow(length(aValsBasic))
  plot(ResultList$proPsRep, ResultList$DiffKS, xlab = "Mean fitness")
  for (i in 1:length(Cols)){
    blnA <- ResultList$aVals == aValsBasic[i]
    points(ResultList$proPsRep[blnA], ResultList$DiffKS[blnA], col = Cols[i])
  }
  legend("bottomright", legend = aValsBasic, col = Cols, pch = 1, cex = 0.5)
  
  # Plot minimum difference per alpha
  plot(ResultList$aVals[idxMinDiffKS], ResultList$DiffKS[idxMinDiffKS], xlab = "alpha")
  
}

PlotResults(ResultList1Full)
PlotResults(ResultList1Fragm)
PlotResults(ResultList2Full)
PlotResults(ResultList2Fragm)

# Function to plot heatmap
PlotHeatmap <- function(ResultList){
  DiffKSMat <- matrix(ResultList$DiffKS, 
                      nrow = length(unique(ResultList$proPsRep)))
  image.plot(x = unique(ResultList$proPsRep),
        y = unique(ResultList$aVals),
        z = DiffKSMat)

}

par(mfrow = c(2, 2))
PlotHeatmap(ResultList2Full_1000G)
PlotHeatmap(ResultList2Fragm_1000G)
heatmap(ResultList2Full_1000G$)
imagehist(FreqFull_1000G)
hist(FreqFragm_1000G)
max(FreqFull_1000G)
max(FreqFragm_1000G)

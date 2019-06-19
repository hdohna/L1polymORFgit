##############################################
#
# General description:
#
#   The following function samples plausible true L1 widths based on simulation
#   results where true and estimated L1 widths were known
#   

# Input:
#
#     SimL1widthTrue: vector of simulated true L1 width values
#     SimL1widthEst:  vector of simulated estimated L1 width values (should
#                     match in length to SimL1widthTrue)
#     SimL1Freq:      vector of population frequency of simulated L1
#     EstL1width:     vector of estimated L1 width (observations)
#     ObsL1Freq:      vector of population frequency of observed L1
#     L1widthBreaks:  vector defining L1 width bins

# Output:
#   
#     SampleMat: rows = true L1 width, column = observed  L1,
#                entries = probability of true L1 width given observed L1


# Comments:
#   
#    Requires package MASS. This function is used in CollapseClosePos_idx

##############################################

SampleTrueL1Width <- function(SimL1widthTrue, 
                              SimL1widthEst,
                              SimL1Freq,
                              EstL1width,
                              ObsL1Freq,
                              L1widthBreaks = seq(0, 6500, 500)
){
  
  # Cut true and estimated L1 width into bins
  SimL1widthTrueCut <- cut(SimL1widthTrue, breaks = L1widthBreaks)
  SimL1widthEstCut  <- cut(SimL1widthEst, breaks  = L1widthBreaks)
  EstL1widthEstCut  <- cut(EstL1width, breaks  = L1widthBreaks)
  
  # Create boolean indicator for estimated and true L1 width to be in the same 
  # bin
  blnSameBin <- SimL1widthTrueCut == SimL1widthEstCut
  
  # Fit logistic regression for probability of estimated length in same bin as
  # true length against estimated frequency
  LogReg <- glm(blnSameBin ~ SimL1Freq, family = binomial)
  summary(LogReg)
  a <- LogReg$coefficients[1]
  b <- LogReg$coefficients[2]
  ExpPart  <- exp(LogReg$coefficients[1] + LogReg$coefficients[2] * ObsL1Freq)
  ProbSame <- ExpPart / (1 + ExpPart)
  plot(ObsL1Freq, ProbSame, col = rgb(0, 0, 0, 0.1), pch = 16)
  
  # Get combinations of L1 width classes among estimated 
  L1combos <- table(SimL1widthEstCut, SimL1widthTrueCut)
  diag(L1combos) <- 0
  
  # Matrix of true L1 width given estimated width
  L1TrueGivenEst <- L1combos / rowSums(L1combos)
  
  # Index matching vector of estimated L1 length to matrix L1TrueGivenEst
  idxRow <- match(EstL1widthEstCut, rownames(L1TrueGivenEst))
  
  # Create sample matrix: rows = true L1 width, column = observed estimated
  # L1 width, entries = probability of true L1 width given estimated
  SampleMat <- sapply(1:length(EstL1widthEstCut), function(x){
    Prob <- L1TrueGivenEst[idxRow[x], ] * (1 - ProbSame[x])
    Prob[idxRow[x]] <- ProbSame[x]
    Prob
  })
  
  # Plot observed vs expected L1 width counts to check whether expected 
  # counts are close to true counts
  plot(rowSums(SampleMat), table(SimL1widthTrueCut), 
       ylab = "True L1 width counts")
  lines(c(0, 1000), c(0, 1000))
  
  # Plot true frequency, estimated frequency, and reconstructed frequency
  # per length bin
  ObsL1Freq_Est    <- aggregate(ObsL1Freq ~ EstL1widthEstCut,     FUN = mean)
  SimL1Freq_True   <- aggregate(SimL1Freq ~ SimL1widthTrueCut, FUN = mean)
  SampledTrueL1idx <- apply(SampleMat, 2, function(x) sample(nrow(SampleMat), 1, prob = x))
  L1MidPts <- 0.5 * (L1widthBreaks[-1] + L1widthBreaks[-length(L1widthBreaks)])
  L1MidPts[SampledTrueL1idx]
}

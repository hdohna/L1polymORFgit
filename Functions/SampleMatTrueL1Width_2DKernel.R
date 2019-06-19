##############################################
#
# General description:
#
#   The following function calculates the probability matrix for a true L1 
#   width (bin) given the estimated L1 width (bin) using a smoothing kernel
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
                              NBins = 12,
                              MaxL1width = 6019
){
  
  # Estimate joint surface for true and estimated L1 width
  PTrueEstKD <- kde2d(x = SimL1widthEst, 
                      y = SimL1widthTrue,
                      n = NBins,
                      lims = c(1, MaxL1width, 1, MaxL1width))
  image(PTrueEstKD)
  
  # Get indices to map L1 width into bins
  idxSimL1widthEst <- sapply(SimL1widthEst, function(x) {
    which.min(abs(x - PTrueEstKD$x))
  })
  idxSimL1widthTrue <- sapply(SimL1widthTrue, function(x) {
    which.min(abs(x - PTrueEstKD$y))
  })
  idxEstL1width <- sapply(EstL1width, function(x) {
    which.min(abs(x - PTrueEstKD$x))
  })

  # Create boolean indicator for estimated and true L1 width to be in the same 
  # bin
  blnSameBin <- idxSimL1widthEst == idxSimL1widthTrue
  
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
  L1combos <- PTrueEstKD$z
  diag(L1combos) <- 0
  
  # Matrix of true L1 width given estimated width
  L1TrueGivenEst <- L1combos / rowSums(L1combos)
  
  # Create sample matrix: rows = true L1 width, column = observed estimated
  # L1 width, entries = probability of true L1 width given estimated
  SampleMat <- sapply(1:length(EstL1widthEstCut), function(x){
    Prob <- L1TrueGivenEst[idxEstL1width[x], ] * (1 - ProbSame[x])
    Prob[idxEstL1width[x]] <- ProbSame[x]
    Prob
  })
  
  # Plot observed vs expected L1 width counts to check whether expected 
  # counts are close to true counts
  SimL1widthTrueCounts <- rep(0, nrow(SampleMat))
  SimL1widthTrueTable  <- table(idxSimL1widthTrue)
  SimL1widthTrueCounts[as.numeric(names(SimL1widthTrueTable))] <- 
    SimL1widthTrueTable
  plot(rowSums(SampleMat), SimL1widthTrueCounts, ylab = "True L1 width counts")
  lines(c(0, 200), c(0, 200))
  
  # Plot true frequency, estimated frequency, and reconstructed frequency
  # per length bin
  ObsL1Freq_Est    <- aggregate(ObsL1Freq ~ idxEstL1width,     FUN = mean)
  SimL1Freq_True   <- aggregate(SimL1Freq ~ idxSimL1widthTrue, FUN = mean)
  SampledTrueL1idx <- apply(SampleMat, 2, function(x) sample(NBins, 1, prob = x))
  ObsL1Freq_Est2   <- aggregate(ObsL1Freq ~ SampledTrueL1idx, FUN = mean)
  
  DiffEstTrue     <- ObsL1Freq_Est$ObsL1Freq - SimL1Freq_True$SimL1Freq
  DiffSampledTrue <- ObsL1Freq_Est2$ObsL1Freq - SimL1Freq_True$SimL1Freq
  plot(PTrueEstKD$x, DiffEstTrue,       col = "blue")
  points(PTrueEstKD$x, DiffSampledTrue, col = "red")
  lines(c(0, 10^4), c(0, 0), lty = 2)
  mean(abs(DiffEstTrue))
  mean(abs(DiffSampledTrue))
  
  SampleMat
}

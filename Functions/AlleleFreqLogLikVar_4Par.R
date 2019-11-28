##############################################
#
# General description:
#
#   The following function calculates the probability log-likelihood of sampled
#   frequencies of alleles under selection (equation after Boissinot at al. 2006 PNAS)
#   accounting for variable L1 length estimation

# Input:
#
#     FreqLengthDF:  data.frame containing
#     L1MidPts: vector of counts midpoints of L1 length classes
#     L1TrueGivenEstList: 
#     Predict: matrix (n x 2) of predictor variable values
#     a:      intercept
#     b:      slope for first predictor
#     c:      slope for second predictor
#     N:      population size
#     DetectProb: probability to detect and insertion

# Comment:
#     This function requires the function AlleleFreqTime and AlleleFreqSample

##############################################

AlleleFreqLogLikVar_4Par <- function(FreqLengthDF, 
                                     L1TrueGivenEstList, 
                                     L1MidPts,
                                     a, 
                                     b, 
                                     c, 
                                     d, 
                                     SD = NULL, 
                                     N, 
                                  SampleSize = 2504,
                                  blnIns, LogRegCoeff,
                                  DetectProb = rep(1, length(Freqs)),
                                  verbose = T, showInfIndices = F,
                                  LowerS = -1,
                                  UpperS = 1){
  # Determine which length classes are full-length
  blnClassFull <- L1MidPts >= LengthFull
  
  # Calculate selection coefficients
  sVals = a + b*L1MidPts + c * blnClassFull
  
  # Calculate integration constants
  IntConsts <- sapply(sVals, function(z){
    AlleleFreqIntConst(s = z, N = N, SampleSize = SampleSize,
                       DetectProb = DetectProb,
                       LogRegCoeff = LogRegCoeff, blnIns = T)
  })
  
  
  LProbs <- sapply(1:nrow(FreqLengthDF), function(i){
    
    # Get index of current length
    idxLength <- FreqLengthDF$L1widthIdx_mean[i]
    
    # Get length values associated with current length and their characteristics
    LengthTrueEst <- L1TrueGivenEstList[[idxLength]] 
    blnFullOther  <- LengthTrueEst$L1Length >= LengthFull
    sValsOther    <- a + b * LengthTrueEst$L1Length + c * blnFullOther
    
    # Calculate the probability of observed length at oberved frequency
    ProbObs <- FreqLengthDF$ProbSameLength_mean[i]  * 
      AlleleFreqSampleVarDiscrete1(k = FreqLengthDF$Freq_mean[i], 
                                   s = sVals[idxLength], 
                                   N = N, IntConst = IntConsts[idxLength], 
                                   SampleSize = SampleSize, DetectProb = DetectProb,
                                   LogRegCoeff = LogRegCoeff, blnIns = T) +
      
      (1 - FreqLengthDF$ProbSameLength_mean[i])  * 
      sum(sapply(1:nrow(LengthTrueEst), function(j){
        LengthTrueEst$Probs[j] *
          AlleleFreqSampleVarDiscrete1(k = FreqLengthDF$Freq_mean[i], 
                                       s = sValsOther[j], 
                                       N = N, IntConst = IntConsts[LengthTrueEst$idxNonZero[j]], 
                                       SampleSize = MEInsSamplesize, DetectProb = DetectProb,
                                       LogRegCoeff = LogRegCoeff, blnIns = T)
      }))
    
    FreqLengthDF$N_sum[i] * log(ProbObs)
    
  })
  if (verbose){
    cat("parameter values: a =", a, "b =", b, "c = ", c, "d = ", d,
        "log-likelihood =", sum(LProbs), "\n")
  }
  
  sum(LProbs)
}



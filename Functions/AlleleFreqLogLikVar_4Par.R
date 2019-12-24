##############################################
#
# General description:
#
#   The following function calculates the probability log-likelihood of sampled
#   frequencies of alleles under selection (equation after Boissinot at al. 2006 PNAS)
#   accounting for variable L1 length estimation

# Input:
#
#     FreqLengthDF:  data.frame containing the following columns
#                    - L1widthIdx_mean: index of current length within length
#                         class vector
#                    - Freq_mean: population frequency
#                    - ProbSameLength_mean: probability that estimated and
#                         true length are the same
#     L1MidPts:      vector of counts midpoints of L1 length classes
#     L1TrueGivenEstList: a list that gives for each estimated length class a 
#                    data.frame that contains the following columns:
#                    - Probs: probabilities of true length classes
#                    - idxNonZero: indices of true length classes
#                    true length classes with non-zero probability and their
#                    probabilities
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
                                     LengthFull = 6000,
                                  SampleSize = 2504,
                                  blnIns, LogRegCoeff,
                                  DetectProb = rep(1, length(Freqs)),
                                  verbose = T, showInfIndices = F,
                                  LowerS = -1,
                                  UpperS = 1,
                                  TaInteraction = c("Intercept", "L1Width", "FullLength")[1]){
  # Determine which length classes are full-length
  blnClassFull <- L1MidPts >= LengthFull
  
  # Calculate selection coefficients
  sVals = a + b*L1MidPts + c * blnClassFull
  sValsTa = switch(TaInteraction, 
                   Intercept = a + d + b * L1MidPts + c * blnClassFull,
                   L1Width = a + (b + d) * L1MidPts + c * blnClassFull,
                   FullLength = a + b * L1MidPts + (c + d) * blnClassFull)
                   
  # Calculate integration constants
  IntConsts <- sapply(sVals, function(z){
    AlleleFreqIntConst(s = z, N = N, SampleSize = SampleSize,
                       DetectProb = DetectProb,
                       LogRegCoeff = LogRegCoeff, blnIns = T)
  })
  
  if (d == 0) { # No difference between Ta and non-Ta LINE-1
    LProbs <- sapply(1:nrow(FreqLengthDF), function(i){
      
      # Get index of current length
      idxLength <- FreqLengthDF$L1widthIdx_mean[i]
      
      # Get length values associated with current length and their characteristics
      LengthTrueEst <- L1TrueGivenEstList[[idxLength]] 
      blnFullOther  <- LengthTrueEst$L1Length >= LengthFull
      sValsOther    <- a + b * LengthTrueEst$L1Length + c * blnFullOther
      
      # Calculate the probability of observed length at oberved frequency
      ProbObs <- 
        
        # Likelihood given that estimated length == true length
        FreqLengthDF$ProbSameLength_mean[i]  * 
        AlleleFreqSampleVarDiscrete1(
          k = FreqLengthDF$Freq_mean[i], 
          s = sVals[idxLength], 
          N = N, IntConst = IntConsts[idxLength], 
          SampleSize = SampleSize, DetectProb = DetectProb,
          LogRegCoeff = LogRegCoeff, blnIns = T) +
        
        # Likelihood given that estimated length != true length
        (1 - FreqLengthDF$ProbSameLength_mean[i])  * 
        sum(sapply(1:nrow(LengthTrueEst), function(j){ # sum over all possible length values
          LengthTrueEst$Probs[j] *
            AlleleFreqSampleVarDiscrete1(
              k = FreqLengthDF$Freq_mean[i], 
              s = sValsOther[j], 
              N = N, IntConst = IntConsts[LengthTrueEst$idxNonZero[j]], 
              SampleSize = SampleSize, DetectProb = DetectProb,
              LogRegCoeff = LogRegCoeff, blnIns = T)
        }))
      
      FreqLengthDF$N_sum[i] * log(ProbObs)
      
    })
    
  } else { # Difference between Ta and non-Ta LINE-1
    IntConstsTa <- sapply(sValsTa, function(z){
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
#      sValsOtherTa  <- a + b * LengthTrueEst$L1Length + (c + d)* blnFullOther
      sValsOtherTa = switch(TaInteraction, 
                       Intercept = a + d + b * LengthTrueEst$L1Length + c * blnFullOther,
                       L1Width = a + (b + d) * LengthTrueEst$L1Length + c * blnFullOther,
                       FullLength = a + b * LengthTrueEst$L1Length + (c + d) * blnFullOther)
      
      # Calculate the probability of observed length at oberved frequency for non-Ta L1
      ProbObs <- FreqLengthDF$ProbSameLength_mean[i]  * 
        AlleleFreqSampleVarDiscrete1(k = FreqLengthDF$Freq_mean[i], 
                                     s = sVals[idxLength], 
                                     N = N, IntConst = IntConsts[idxLength], 
                                     SampleSize = SampleSize, DetectProb = DetectProb,
                                     LogRegCoeff = LogRegCoeff, blnIns = T) +
        
        (1 - FreqLengthDF$ProbSameLength_mean[i])  * 
        sum(sapply(1:nrow(LengthTrueEst), function(j){
          LengthTrueEst$Probs[j] *
            AlleleFreqSampleVarDiscrete1(
              k = FreqLengthDF$Freq_mean[i], 
              s = sValsOther[j], 
              N = N, IntConst = IntConsts[LengthTrueEst$idxNonZero[j]], 
              SampleSize = SampleSize, DetectProb = DetectProb,
              LogRegCoeff = LogRegCoeff, blnIns = T)
        }))
      
      # Calculate the probability of observed length at oberved frequency for non-Ta L1
      ProbObsTa <- FreqLengthDF$ProbSameLength_mean[i]  * 
        AlleleFreqSampleVarDiscrete1(k = FreqLengthDF$Freq_mean[i], 
                                     s = sValsTa[idxLength], 
                                     N = N, IntConst = IntConstsTa[idxLength], 
                                     SampleSize = SampleSize, DetectProb = DetectProb,
                                     LogRegCoeff = LogRegCoeff, blnIns = T) +
        
        (1 - FreqLengthDF$ProbSameLength_mean[i])  * 
        sum(sapply(1:nrow(LengthTrueEst), function(j){
          LengthTrueEst$Probs[j] *
            AlleleFreqSampleVarDiscrete1(k = FreqLengthDF$Freq_mean[i], 
                                         s = sValsOtherTa[j], 
                                         N = N, IntConst = IntConsts[LengthTrueEst$idxNonZero[j]], 
                                         SampleSize = SampleSize, DetectProb = DetectProb,
                                         LogRegCoeff = LogRegCoeff, blnIns = T)
        }))
      (FreqLengthDF$N_sum[i] - FreqLengthDF$blnTa_sum[i]) * log(ProbObs) + 
        FreqLengthDF$blnTa_sum[i] * log(ProbObsTa)
      
    })
    
  }
  if (verbose){
    cat("parameter values: a =", a, "b =", b, "c = ", c, "d = ", d,
        "log-likelihood =", sum(LProbs), "\n")
  }
  
  sum(LProbs)
}



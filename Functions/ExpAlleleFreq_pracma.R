##############################################
#
# General description:
#
#   The following function calculates the expected frequency of an allele
#   under selection (equation after Boissinot at al. 2006 PNAS)
#   (calculation is conditional on inclusion in study)

# Input:
#
#     s: selection coefficient
#     N: population size
#     SampleSize: sample size

# Comment:
#     This function requires the function AlleleFreqTime

##############################################

ExpAlleleFreq_pracma <- function(s, N, SampleSize = 2504,
                                  blnIns = T, LogRegCoeff,
                                  DetectProb = rep(1, length(Freqs))){
    
  CountVals  <- 1:(SampleSize - 1)
  CountProbs <- sapply(CountVals, function(i){
    exp(AlleleFreqSample_pracma(k = i, s = s, N = N, SampleSize = SampleSize, blnIns = blnIns,
                         LogRegCoeff = LogRegCoeff,
                         DetectProb = DetectProb))
    
  })
  (CountVals %*% CountProbs) / SampleSize / sum(CountProbs)
}





##############################################
#
# General description:
#
#   The following function calculates the probability log-likelihood of sampled
#   frequencies of alleles under selection (equation after Boissinot at al. 2006 PNAS)

# Input:
#
#     Freqs:  vector of observed allele frequencies (integer)
#     Counts: vector of counts (number of alleles with the frequencies 
#             specified in Freqs)
#     Predict: matrix (n x 2) of predictor variable values
#     a:      intercept
#     b:      slope for first predictor
#     c:      slope for second predictor
#     N:      population size

# Comment:
#     This function requires the function AlleleFreqTime and AlleleFreqSample

##############################################

AlleleFreqLogLik_4Par <- function(Freqs, Counts, Predict, a, b, c, d, N, SampleSize = 2504,
                                 MinFactor = 2, blnUseFPrime = T){
  if ((length(Freqs) != length(Counts)) | (length(Freqs) != nrow(Predict)) |
      (nrow(Predict) != length(Counts)) ){
    stop("Freqs, Counts and Predict vector have to have the same length\n")
  }
  LogLikVals <- sapply(1:length(Freqs), function(i){
    Counts[i] * AlleleFreqSample(Freqs[i], a + b * Predict[i, 1] +
                                       c * Predict[i, 2] + d * Predict[i, 3], N, 
                                 SampleSize = SampleSize, blnUseFPrime = blnUseFPrime)
  })
  blnInf <- is.infinite(LogLikVals)
  LogLikVals[blnInf] <- MinFactor * min(LogLikVals[!blnInf])
  sum(LogLikVals)
}



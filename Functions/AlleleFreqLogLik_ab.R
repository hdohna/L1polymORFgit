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
#     Predict: vector of predictor variable values
#     s:      selection coefficient
#     N:      population size

# Comment:
#     This function requires the function AlleleFreqTime and AlleleFreqSample

##############################################

AlleleFreqLogLik_ab <- function(Freqs, Counts, Predict, a, b, N, SampleSize = 2504){
  if ((length(Freqs) != length(Counts)) | (length(Freqs) != length(Predict)) |
      (length(Predict) != length(Counts)) ){
    stop("Freqs, Counts and Predict vector have to have the same length\n")
  }
  sum(sapply(1:length(Freqs), function(i){
    Counts[i] * log(AlleleFreqSample(Freqs[i], a + b * Predict[i], N, 
                                 SampleSize = SampleSize))
  }))
  
}



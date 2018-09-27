##############################################
#
# General description:
#
#   The following function calculates the derivative of log-likelihood 
#   of sampled
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

AlleleFreqLogLik_abc_Deriv <- function(Freqs, Counts, Predict, a, b, c, N, SampleSize = 2504,
                                 MinFactor = 2){
  if ((length(Freqs) != length(Counts)) | (length(Freqs) != nrow(Predict)) |
      (nrow(Predict) != length(Counts)) ){
    stop("Freqs, Counts and Predict vector have to have the same length\n")
  }
  LogLikDer <- sapply(1:length(Freqs), function(i){
    Counts[i] * AlleleFreqSample_Deriv(Freqs[i], a + b * Predict[i, 1] +
                                       c* Predict[i, 2], N, 
                                 SampleSize = SampleSize)
  })
  c(a = sum(LogLikDer), 
    b = Predict[i, 1] %*% LogLikDer, 
    c = Predict[i, 2] %*% LogLikDer)
}



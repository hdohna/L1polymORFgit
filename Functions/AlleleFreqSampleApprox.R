##############################################
#
# General description:
#
#   The following function calculates the log probability of the count value
#   for an allele under selection (equation after Boissinot at al. 2006 PNAS)

# Input:
#
#     k: number of samples carrying the allele
#     s: selection coefficient
#     N: population size
#     SampleSize: sample size

# Comment:
#     This function requires the function AlleleFreqTime

##############################################

AlleleFreqSampleApprox <- function(k, s, N, SampleSize = 2504, blnUseFPrime = T,
                                   NrPts = 100,
                                   MaxW = 0.1){
    
  # The lines below use F' in Boissinot's paper
  if (blnUseFPrime){
    
    # Calculate integration constant
    IntConst <- integrate(function(x) AlleleFreqTime(x, s, N), 
                          0, 1)$value - 
      integrate(function(x) (1 - x)^SampleSize * AlleleFreqTime(x, s, N), 
                0, 1)$value
    
    
    # Calculate probability of obtaining k alleles in a sample of size 
    # SampleSize
    xVals <- seq(0, 1, 1/NrPts)
    xVals <- xVals[-length(xVals)] + 0.5/NrPts
    LChoose <- lchoose(SampleSize, k)
    LogVals <- sapply(xVals, function(x){
      log(1 - (1 - x)^SampleSize) + log(AlleleFreqTime(x, s, N)) + k*log(x) +
        (SampleSize - k) * log(1 - x)
    }) + LChoose
    sum(LogVals) - (NrPts - 1) * (MaxW *max(LogVals) + (1- MaxW) *min(LogVals)) - log(IntConst)
  } else {
    
    # Calculate integration constant
    IntConst <- integrate(function(x) AlleleFreqTime(x, s, N), 0, 1)$value
    
    # Calculate probability of obtaining k alleles in a sample of size 
    # SampleSize
    lchoose(SampleSize, k) +
      log(integrate(function(x) AlleleFreqTime(x, s, N) * x^k * 
                      (1 - x)^(SampleSize - k), 0, 1)$value) - log(IntConst)
  }
  
  
}



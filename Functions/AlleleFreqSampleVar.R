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
#     This function requires the package pracma

##############################################

AlleleFreqSampleVar <- function(k, m, SD, N = 10^4, SampleSize = 2*2504, LowerS = -1,
                                UpperS = 1){
    
  # Calculate integration constant
  IntConst <- integral2(fun = function(x, y){
     dnorm(y, m, SD) * AlleleFreqTime(x, y, N) * 
       (1 - (1 - x)^SampleSize - x^SampleSize)
   }, xmin = 0, xmax = 1, ymin = LowerS, ymax = UpperS)$Q
  
    
    # Calculate probability of obtaining k alleles in a sample of size 
    # SampleSize
    log(integral2(fun = function(x, y){
        dnorm(y, m, SD) * AlleleFreqTime(x, y, N) * 
        (1 - (1 - x)^SampleSize - x^SampleSize) *
        dbinom(k, SampleSize, x)
      }, xmin = 0, xmax = 1, ymin = LowerS, ymax = UpperS)$Q
    ) - log(IntConst)
}



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

AlleleFreqSampleVar <- function(k, m, SD, N = 10^4, SampleSize = 2*2504, 
                                LogRegCoeff,
                                DetectProb = 1,
                                blnIns = T,
                                LowerS = -1,
                                UpperS = 1){
    
  ProbRef <- function(x) {
    ExpCoeff <- exp(LogRegCoeff[1] + LogRegCoeff[2] * x)
    ExpCoeff / (1 + ExpCoeff)
  }
  
  # The lines below are for insertions relative to the reference
  if (blnIns){
    
    # Calculate integration constant
    IntConst <- integral2(fun = function(x, y){
      dnorm(y, m, SD) * AlleleFreqTime(x, y, N) * 
        (1 - (1 - x)^SampleSize) * (1 - ProbRef(x))
    }, xmin = 0, xmax = 1, ymin = LowerS, ymax = UpperS)$Q
    
    
    # Calculate probability of obtaining k alleles in a sample of size 
    # SampleSize
    log(integral2(fun = function(x, y){
      dnorm(y, m, SD) * AlleleFreqTime(x, y, N) * 
        (1 - (1 - x)^SampleSize - x^SampleSize) *
        dbinom(k, SampleSize, DetectProb * x) * (1 - ProbRef(x))
    }, xmin = 0, xmax = 1, ymin = LowerS, ymax = UpperS)$Q
    ) - log(IntConst)
    
  } else {
    
    # Calculate integration constant
    IntConst <- integral2(fun = function(x, y){
      dnorm(y, m, SD) * AlleleFreqTime(x, y, N) * 
        (1 - x^SampleSize) * ProbRef(x)
    }, xmin = 0, xmax = 1, ymin = LowerS, ymax = UpperS)$Q
    
    
    # Calculate probability of obtaining k alleles in a sample of size 
    # SampleSize
    log(integral2(fun = function(x, y){
      dnorm(y, m, SD) * AlleleFreqTime(x, y, N) * 
        (1 - x^SampleSize) *
        dbinom(k, SampleSize, DetectProb * (1 - x)) * ProbRef(x)
    }, xmin = 0, xmax = 1, ymin = LowerS, ymax = UpperS)$Q
    ) - log(IntConst)
    
  }
  
 }



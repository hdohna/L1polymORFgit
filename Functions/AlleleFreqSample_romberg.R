##############################################
#
# General description:
#
#   The following function calculates the log probability of the count value
#   for an allele under selection (equation after Boissinot at al. 2006 PNAS)
#   using Romberg integration
# Input:
#
#     k: number of samples carrying the allele
#     s: selection coefficient
#     N: population size
#     SampleSize: sample size
#     DetectProb: probability to detect insertion
#     LogRegCoeff: vector (of length 2) of regression coefficients indicationg
#        how the probability to be on the reference genome depends on allele
#        frequency
#     blnIns: boolean indicator whether analyzing L1 insertions or L1 deletions 

# Comment:
#     This function requires the function AlleleFreqTime and package pracma

##############################################

AlleleFreqSample_romberg <- function(k, s, N, SampleSize = 2504, DetectProb = 1,
                             LogRegCoeff, blnIns = T){
    
  ProbRef <- function(x) {
    ExpCoeff <- exp(LogRegCoeff[1] + LogRegCoeff[2] * x)
    ExpCoeff / (1 + ExpCoeff)
  }
  if (blnIns){  
    AFTime <- function(x, s, N, SampleSize, DetectProb) {
      AlleleFreqTime(x, s, N) * (1 - ProbRef(x))
    }
    AFTimeSample <- function(x, s, N, SampleSize, DetectProb) {
      AlleleFreqTime(x, s, N) * (1 - DetectProb * x)^SampleSize * (1 - ProbRef(x))
    }
    AFTimeBinom <- function(x, s, N, k, SampleSize, DetectProb) {
      AlleleFreqTime(x, s, N) * (1 - ProbRef(x)) * 
        dbinom(k, SampleSize, DetectProb * x)
    }
    AFTimeSampleBinom <- function(x, s, N, k, SampleSize, DetectProb) {
      AlleleFreqTime(x, s, N) * (1 - DetectProb * x)^SampleSize * (1 - ProbRef(x)) * 
        dbinom(k, SampleSize, DetectProb * x)
    }
    
  } else {
    AFTime <- function(x, s, N, SampleSize, DetectProb) {
      AlleleFreqTime(x, s, N) * ProbRef(x)
    }
    AFTimeSample <- function(x, s, N, SampleSize, DetectProb) {
      AlleleFreqTime(x, s, N) * (1 - (1-x)*DetectProb)^SampleSize * ProbRef(x)
    }
    AFTimeBinom <- function(x, s, N, k, SampleSize, DetectProb) {
      AlleleFreqTime(x, s, N) * ProbRef(x) * 
        dbinom(k, SampleSize, DetectProb * (1 - x))
    }
    AFTimeSampleBinom <- function(x, s, N, k, SampleSize, DetectProb) {
      AlleleFreqTime(x, s, N) * (1 - (1-x)*DetectProb)^SampleSize * ProbRef(x) * 
        dbinom(k, SampleSize, DetectProb * (1 - x))
    }
    
  }

     # Calculate integration constant
    IntConst <- romberg(AFTime,  0, 1, s = s, N = N, SampleSize = SampleSize, 
                        DetectProb = DetectProb)$value - 
    romberg(AFTimeSample,  0, 1, s = s, N = N, SampleSize = SampleSize, 
            DetectProb = DetectProb)$value
    
    # Calculate probability of obtaining k alleles in a sample of size 
    # SampleSize
    log(
      romberg(AFTimeBinom, 0, 1, s = s, N = N, k = k, SampleSize = SampleSize, 
              DetectProb = DetectProb)$value - 
      romberg(AFTimeBinom, 0, 1, s = s, N = N, k = k, SampleSize = SampleSize, 
              DetectProb = DetectProb)$value
      ) - log(IntConst)
  
} 


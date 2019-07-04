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
#     DetectProb: probability to detect insertion
#     LogRegCoeff: vector (of length 2) of regression coefficients indicationg
#        how the probability to be on the reference genome depends on allele
#        frequency
#     blnIns: boolean indicator whether analyzing L1 insertions or L1 deletions 

# Comment:
#     This function requires the function AlleleFreqTime

##############################################

AlleleFreqSample_simplified2 <- function(k, s, N, SampleSize = 2504, DetectProb = 1,
                             LogRegCoeff, blnIns = T){
    
  # The lines below are for insertions relative to the reference
  if (blnIns){
    
    # # Calculate integration constant
    # IntConst <- integrate(function(x) {
    #      AlleleFreqTime(x, s, N) * (1 - (1 - x)^SampleSize) * (1 - ProbRef(x))
    #   },  0, 1)$value
    # 
    # # Calculate probability of obtaining k alleles in a sample of size 
    # # SampleSize
    # log(integrate(function(x) {
    #   AlleleFreqTime(x, s, N) * (1 - (1 - x)^SampleSize) * (1 - ProbRef(x)) * 
    #                 dbinom(k, SampleSize, DetectProb * x)
    #   }, 0, 1)$value) -
    # log(IntConst)
  
    # Calculate integration constant
    IntConst <- integrate(function(x) {
        AlleleFreqTime(x, s, N) 
      },  0, 1)$value - 
    integrate(function(x) {
        AlleleFreqTime(x, s, N) * (1 - DetectProb * x)^SampleSize 
      },  0, 1)$value
    
    # Calculate probability of obtaining k alleles in a sample of size 
    # SampleSize
    log(
      integrate(function(x) {
          AlleleFreqTime(x, s, N) * 
          dbinom(k, SampleSize, DetectProb * x)
        }, 0, 1)$value - 
      integrate(function(x) {
          AlleleFreqTime(x, s, N) * (1 - DetectProb * x)^SampleSize * 
            dbinom(k, SampleSize, DetectProb * x)
        }, 0, 1)$value
      ) -
      log(IntConst)
    # The lines below are for deletions relative to the reference
  } else {
    
    # # Calculate integration constant
    # IntConst <- integrate(function(x) {
    #   AlleleFreqTime(x, s, N) * (1 - x^SampleSize) * ProbRef(x)
    # },  0, 1)$value
    # 
    # # Calculate probability of obtaining k alleles in a sample of size 
    # # SampleSize
    # log(integrate(function(x) {
    #     AlleleFreqTime(x, s, N) * (1 - x^SampleSize) * ProbRef(x) * 
    #     dbinom(k, SampleSize, DetectProb * (1 - x))
    #   }, 0, 1)$value) - log(IntConst)

        # Calculate integration constant
    # IntConst <- integrate(function(x) {
    #      AlleleFreqTime(x, s, N) * ProbRef(x)
    #   },  0, 1)$value - 
    # integrate(function(x) {
    #     AlleleFreqTime(x, s, N) * (1 - (1-x)*DetectProb)^SampleSize * ProbRef(x)
    #   },  0, 1)$value
    # 
    # # Calculate probability of obtaining k alleles in a sample of size 
    # # SampleSize
    # log(
    #   integrate(function(x) {
    #       AlleleFreqTime(x, s, N) * ProbRef(x) * 
    #       dbinom(k, SampleSize, DetectProb * (1 - x))
    #     }, 0, 1)$value - 
    #   integrate(function(x) {
    #     AlleleFreqTime(x, s, N) * (1 - (1-x)*DetectProb)^SampleSize * ProbRef(x) * 
    #       dbinom(k, SampleSize, DetectProb * (1 - x))
    #   }, 0, 1)$value
    # ) - log(IntConst)
  }
  
  
}



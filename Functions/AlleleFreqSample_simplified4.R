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

AlleleFreqSample_simplified4 <- function(k, s, N, SampleSize = 2504, DetectProb = 1,
                             LogRegCoeff, blnIns = T){
    

    # Calculate integration constant
    IntConst <- integral(function(x) {
        AlleleFreqTime(x, s, N) 
      },  0, 1)
    
    # Calculate probability of obtaining k alleles in a sample of size 
    # SampleSize
    #log(
    integral(function(x) {
          AlleleFreqTime(x, s, N) * 
          dbinom(k, SampleSize, x)
        }, 0, 1) /
      #) -
      IntConst
  }



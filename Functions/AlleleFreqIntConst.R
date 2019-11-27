##############################################
#
# General description:
#
#   The following function calculates the integration constant necessary to 
#   for calculationg the  probability of the count value
#   for an allele under selection (equation after Boissinot at al. 2006 PNAS)

# Input:
#
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

AlleleFreqIntConst <- function(s, N, SampleSize = 2504, DetectProb = 1,
                             LogRegCoeff, blnIns = T){
    
  ProbRef <- function(x) {
    ExpCoeff <- exp(LogRegCoeff[1] + LogRegCoeff[2] * x)
    ExpCoeff / (1 + ExpCoeff)
  }
  # The lines below are for insertions relative to the reference
  if (blnIns){
    
    # Calculate integration constant
    IntConst <- integral(function(x) {
        AlleleFreqTime(x, s, N) * (1 - ProbRef(x))
      },  0, 1) - 
    integral(function(x) {
        AlleleFreqTime(x, s, N) * (1 - DetectProb * x)^SampleSize * (1 - ProbRef(x))
      },  0, 1)
    
   } else {
    
    # Calculate integration constant
    IntConst <- integral(function(x) {
         AlleleFreqTime(x, s, N) * ProbRef(x)
      },  0, 1) - 
    integral(function(x) {
        AlleleFreqTime(x, s, N) * (1 - (1-x)*DetectProb)^SampleSize * ProbRef(x)
      },  0, 1)
    
   }
  IntConst
  
}



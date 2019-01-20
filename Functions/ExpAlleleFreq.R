##############################################
#
# General description:
#
#   The following function calculates the expected frequency of an allele
#   under selection (equation after Boissinot at al. 2006 PNAS)
#   (calculation is conditional on inclusion in study)

# Input:
#
#     s: selection coefficient
#     N: population size
#     SampleSize: sample size

# Comment:
#     This function requires the function AlleleFreqTime

##############################################

ExpAlleleFreq <- function(s, N = 10^4, SampleSize = 2*2504, blnIns = T){
    
    # # Calculate integration constant
    # integrate(function(x) x * SampleSize / 
    #             (1 - (1 - x)^SampleSize - x^SampleSize) * AlleleFreqTime(x, s, N), 0, 1)$value /
    # integrate(function(x) AlleleFreqTime(x, s, N), 0, 1)$value

  ProbRef <- function(x) {
    ExpCoeff <- exp(LogRegCoeff[1] + LogRegCoeff[2] * x)
    ExpCoeff / (1 + ExpCoeff)
  }
  
  if (blnIns) {
    IntConst <- integrate(function(x) (1 - (1 - x)^SampleSize) * ProbRef(x) *
                            AlleleFreqTime(x, s, N), 
                          0, 1)$value
    # Calculate expected frequency
    integrate(function(x) x* (1 - (1 - x)^SampleSize) * ProbRef(x) *
                AlleleFreqTime(x, s, N) * x, 0, 1)$value / IntConst
    
  }
  # Calculate integration constant
}


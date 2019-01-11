##############################################
#
# General description:
#
#   The following function calculates the posterior probability density of the
#   frequency of an allele under selection (equation after Boissinot at al. 2006 PNAS)

# Input:
#
#     y: allele frequency
#     s: selection coefficient
#     N: population size

# Comment:
#     This function requires the function AlleleFreqTime

##############################################

AlleleFreqPostProb <- function(y, s, N = 10^4, blnIns = T){
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
}



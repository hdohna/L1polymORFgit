##############################################
#
# General description:
#
#   The following function calculates the probability that an allele is 
#   becomes fixed at proportion 1 (equation 5.47 in Ewens 2004 Mathematical 
#   Population Genetics)

# Input:
#
#     s: selection coefficient
#     N: population size

# Comment:
#     This function requires the function AlleleFreqTime

##############################################

ProbFix1 <- function(s, N = 10^4){
  if (s != 0){
    (1 - exp(-s)) / (1 - exp(-2*N*s))
  } else {
    exp(-s)/(2*N*exp(-2*N*s))
  }
}



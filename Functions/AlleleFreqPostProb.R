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

AlleleFreqPostProb <- function(y, s, N){
  y * AlleleFreqTime(y, s, N) / 
    integrate(function(x) x * AlleleFreqTime(x, s, N), 0, 1)$value
}



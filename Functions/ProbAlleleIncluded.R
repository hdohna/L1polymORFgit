##############################################
#
# General description:
#
#   The following function calculates the probability that an allele is 
#   included in a study
#   for an allele under selection (equation after Boissinot at al. 2006 PNAS)

# Input:
#
#     s: selection coefficient
#     N: population size
#     SampleSize: sample size

# Comment:
#     This function requires the function AlleleFreqTime

##############################################

ProbAlleleIncluded <- function(s, N = 10^4, SampleSize = 2*2504){
    
    # Calculate integration constant
    integrate(function(x)  (1 - (1 - x)^SampleSize - x^SampleSize) * 
                AlleleFreqTime(x, s, N), 0, 1)$value /
    integrate(function(x) AlleleFreqTime(x, s, N), 0, 1)$value
}



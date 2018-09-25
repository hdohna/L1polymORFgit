##############################################
#
# General description:
#
#   The following function calculates the probability density of the frequency 
#   of an allele under selection (equation after Boissinot at al. 2006 PNAS)

# Input:
#
#     k: number of samples carrying the allele
#     s: selection coefficient
#     N: population size
#     SampleSize: sample size

# Comment:
#     This function requires the function AlleleFreqTime

##############################################

AlleleFreqSample <- function(k, s, N, SampleSize = 2504){
    
   # Calculate integration constant
   IntConst <- integrate(function(x) x * AlleleFreqTime(x, s, N), 0, 1)$value
   
   # Calculate probability of obtaining k alleles in a sample size 
   choose(SampleSize, k) * 
   integrate(function(x) AlleleFreqTime(x, s, N) * x^(k + 1) *
                (1 - x)^(SampleSize - k) , 0, 1)$value /
    IntConst
}



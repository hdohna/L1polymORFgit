##############################################
#
# General description:
#
#   The following function calculates the  log probability of all count values
#   for an allele under selection (equation after Boissinot at al. 2006 PNAS)

# Input:
#
#     s: selection coefficient
#     N: population size
#     SampleSize: sample size

# Comment:
#     This function requires the function AlleleFreqTime

##############################################

AlleleFreqSampleProb <- function(s, N, SampleSize = 2504, MinFactor = 2){
    
   # Calculate integration constant
#   IntConst <- integrate(function(x) x * AlleleFreqTime(x, s, N), 0, 1)$value
   
   # Calculate probability of obtaining k alleles in a sample size
   sapply(0:SampleSize, function(k){
     lchoose(SampleSize, k) + 
       integrate(function(x) log(AlleleFreqTime(x, s = 0, N = 10^4)) + 
                   k * log(x) + (SampleSize - k) * log(1 - x), 0, 1)$value       
   }) #- log(IntConst)
}



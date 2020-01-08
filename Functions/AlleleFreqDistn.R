##############################################
#
# General description:
#
#   The following function calculates the probability density of the
#   frequency of an allele under selection (equation after Boissinot at al. 2006 PNAS)

# Input:
#
#     y: allele frequency
#     s: selection coefficient
#     N: population size

# Comment:
#     This function requires the function AlleleFreqTime

##############################################

AlleleFreqDistn <- function(y, s, N = 10^4){
    
  IntConst <- integral(function(x) AlleleFreqTime(x, s, N), 0, 1)
  
  # Calculate expected frequency
  AlleleFreqTime(y, s, N) / IntConst
    
}



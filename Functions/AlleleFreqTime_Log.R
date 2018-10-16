##############################################
#
# General description:
#
#   The following function calculates the time an allele under selection spends
#   at a particular frequency (equation after Boissinot at al. 2006 PNAS)

# Input:
#
#     y: allele frequency
#     s: selection coefficient
#     N: population size

# Comment:
#     This function is used in the function AlleleFreqSample

##############################################

AlleleFreqTime_Log <- function(y, s, N) {
  if (s != 0){
    -s + (1 + sign(1/(2*N) - y)) / 2 * (log(exp(s) - exp(2*N*s)) + log((exp(2*N*s*y) - 1))) -
             - (1 + sign(y - 1/(2*N))) / 2  * log(exp(s) - 1) + log(exp(2*N*s) - exp(2*N*s*y)) - 
      log((exp(2*N*s) - 1)*N*s*(y - 1)*y)
    
  } else {
    (1 + sign(1/(2*N) - y)) / 2 * log(1 - 2*N) - log(N*(y - 1)*(N + 1)) -
    (1 + sign(y - 1/(2*N))) / 2 * log(N*y*(N + 1))
    
  }
  
}



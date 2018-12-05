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
#     This function is used in the function AlleleFreqSampleVar
#     Requires package pracma

##############################################

AlleleFreqTimeVar <- function(y, s, SD, m, N, DetectProb = 1,
                              blnIns = T) {
  # if (s != 0){
  #   exp(-s)*(((1 + sign(1/(2*N) - y)) / 2 * (exp(s) - exp(2*N*s))*(exp(2*N*s*y) - 1))
  #            - (1 + sign(y - 1/(2*N))) / 2  * ((exp(s) - 1)*(exp(2*N*s) - exp(2*N*s*y)))) / 
  #     ((exp(2*N*s) - 1)*N*s*(y - 1)*y)
  #   
  # } else {
  #   (1 + sign(1/(2*N) - y)) / 2 * (1 - 2*N) / (N*(y - 1)) +
  #   (1 + sign(y - 1/(2*N))) / 2 / (N*y)
  #   
  # }
  NormF <- exp(-(s - m)^2 / SD)
  AFTime <-
    exp(-s)*(((1 + sign(1/(2*N) - y)) / 2 * (exp(s) - exp(2*N*s))*(exp(2*N*s*y) - 1))
             - (1 + sign(y - 1/(2*N))) / 2  * ((exp(s) - 1)*(exp(2*N*s) - exp(2*N*s*y)))) / 
      ((exp(2*N*s) - 1)*N*s*(y - 1)*y)
  AddTerm <- is.na(AFTime) *((1 + sign(1/(2*N) - y)) / 2 * (1 - 2*N) / (N*(y - 1)) +
      (1 + sign(y - 1/(2*N))) / 2 / (N*y))
  AFTime[is.na(AFTime)] <- 0 
  NormF * (AFTime + AddTerm)
}



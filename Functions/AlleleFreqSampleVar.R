##############################################
#
# General description:
#
#   The following function calculates the log probability of the count value
#   for an allele under selection (equation after Boissinot at al. 2006 PNAS)

# Input:
#
#     k: number of samples carrying the allele
#     s: selection coefficient
#     N: population size
#     SampleSize: sample size

# Comment:
#     This function requires the function AlleleFreqTime

##############################################

AlleleFreqSampleVar <- function(k, m, SD, N, SampleSize = 2504, blnUseFPrime = T){
    
  # The lines below use F' in Boissinot's paper
  if (blnUseFPrime){
    
    # Calculate integration constant
    IntConst <- integrate(function(s){
      dnorm(s, m, SD) * integrate(function(x) AlleleFreqTime(x, s, N), 
                0, 1)$value - 
        integrate(function(x) (1 - x)^SampleSize * AlleleFreqTime(x, s, N), 
                  0, 1)$value
      
    }, -1, 1)$value
    
    
    # Calculate probability of obtaining k alleles in a sample of size 
    # SampleSize
    # lchoose(SampleSize, k) +
    #   log(integrate(function(x) AlleleFreqTime(x, s, N) * x^k * 
    #                   (1 - x)^(SampleSize - k), 0, 1)$value -
    #   integrate(function(x) (1 - x)^SampleSize * AlleleFreqTime(x, s, N) * x^k * 
    #   (1 - x)^(SampleSize - k), 0, 1)$value) - log(IntConst)

      log(integrate(function(x) AlleleFreqTime(x, s, N)*dbinom(k, SampleSize, x), 
                    0, 1)$value -
            integrate(function(x) (1 - x)^SampleSize * AlleleFreqTime(x, s, N)*
                        dbinom(k, SampleSize, x), 0, 1)$value) - log(IntConst)
  } else {
    
    # Calculate integration constant
    IntConst <- integrate(function(x) AlleleFreqTime(x, s, N), 0, 1)$value
    
    # Calculate probability of obtaining k alleles in a sample of size 
    # SampleSize
    log(integrate(function(x) AlleleFreqTime(x, s, N) *
                      dbinom(k, SampleSize, x), 0, 1)$value) - log(IntConst)
  }
  
  
}



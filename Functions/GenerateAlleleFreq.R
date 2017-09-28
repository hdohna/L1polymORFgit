##############################################
#
# General description:
#
#   The following function generates allele frequencies using a simple 
#   selection model
#   

# Input:
#
#     Gshape: shape parameter of gamma distribution of fitness effects
#     GSscale: scale parameter of gamma distribution of fitness effects
#     n: population size
#     NrGen: number of generations
#     NrRep: number of alleles generated
#     NrSamples: number of alleles sampled in the end

# Output:
#   
#     SampledFreq: simulated sampled allele frequencies

##############################################

GenerateAlleleFreq <- function(Gshape, GSscale, n = 10^4, NrGen = 10^3,
                               NrRep = 10^4, NrSamples = 10^3){
  cat("Simulating allele frequencies with Gshape =", Gshape, 
      "and GSscale = ", GSscale, " ....")
  
  # Create a vector of selection coefficients
  PCoeff  <- rgamma(NrRep, shape = Gshape, scale = GSscale)
  
  # Initialize vector of allele frequencies (all start with initial frequency
  # 1/2n)
  AlleleFreq <- rep(1/(2*n), length(PCoeff))
  
  # Loop over generations and calculate new allele frequencies
  for (i in 1:NrGen){
    #if (any(is.na(AlleleFreq))) browser()
    SampleProbs <- AlleleFreq * PCoeff / (1 + AlleleFreq * (PCoeff - 1) )
    AlleleFreq  <- rbinom(length(SampleProbs), 2*n, SampleProbs) / (2*n)
    #    AlleleFreq  <- pmin(1 - 1/(2*n), AlleleFreq)
    AlleleFreq  <- pmax(1/(2*n), AlleleFreq)
  }
  cat("done!\n")
  SampledFreq <- sample(AlleleFreq, NrSamples, prob = AlleleFreq)
}

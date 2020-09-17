##############################################
#
# General description:
#
#   The following function simulates a site-frequency spectrum (SFS) based
#   on a selection coefficient and calculates the multinomial probabilities

# Input:
#
#     s: selection coefficient
#     SampleSize: number of alleles sampled
#     N: population size
#     NGen: number of generations to simulate

# Comment:

##############################################

SimulateSFS <- function(s = -0.001, SampleSize = 5000,
                        NMutations = 10^4,
                        N = 10^4, NGen = 5*10^4,
                        PDiffThresh = 10^(-3)){
    

  # Initialize alleles
  Alleles <- 1:(2*N)
  
  
   # Initialize vector of allele frequencies
  PCoeff      <- 1 + s
  AlleleFreq  <- rep(1/(2*N), NMutations)
  MeanFreqOld <- 1
  PDiff       <- 1
  MeanV       <- rep(NA, NGen)
  
  # Loop over generations and calculate new allele frequencies
  cat("Running simulation for", NGen, "generations ... ")
  for(i in 1:NGen){
    
    SampleProbs <- AlleleFreq * PCoeff / (1 + AlleleFreq * (PCoeff - 1) )
    AlleleFreq  <- rbinom(length(SampleProbs), 2*N, SampleProbs)/(2*N)
    MeanFreq    <- mean(AlleleFreq)
    AlleleFreq[AlleleFreq == 0] <- 1/(2*N)
    AlleleFreq[AlleleFreq == 1] <- 1/(2*N)
    # PDiff <- abs(MeanFreq - MeanFreqOld)/MeanFreq
    # MeanFreqOld <- MeanFreq
    MeanV[i] <- MeanFreq

  }
  cat("done!\n")
  hist(AlleleFreq * NMutations, breaks = 0:(2*N), xlim = c(0, 30))
  plot(MeanV, type = "l")
  table(AlleleFreq)
  AlleleFreq * SampleSize
}



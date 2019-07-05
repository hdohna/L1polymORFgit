##############################################
#
# General description:
#
#   The following function simulates a site-frequency spectrum (SFS) based
#   on a selection coefficient and calculates the multinomial probabilities

# Input:
#
#     y: allele frequency
#     s: selection coefficient
#     N: population size

# Comment:
#     This function requires the function AlleleFreqTime

##############################################

SimulateMultinomSFSSelection <- function(SFS, FreqCountV, s, SampleSize = 5000,
                                         N = 10^4, NAdd = 100, 
                                         NGenPerBatch = 10^3, MaxNGen = 10^4,
                                         PDiffThresh = 0.3){
    
  # Check whether SFS has one more element than FreqCountV (the last element
  # in SFS is the number of counts larger than max in SFS)
  if ((length(FreqCountV) + 1) != length(SFS)){
    stop("Input vectors SFS has to have one more element than FreqCountV!")
  }
  
  # Initialize vector of allele frequencies
  PCoeff <- 1 + s
  AlleleFreqStart <- rep(1/(2*N), NAdd)
  AlleleFreq      <- AlleleFreqStart
  NDropOut        <- 0
  MeanFreq        <- 0
  NrAlleles       <- length(AlleleFreq)
  PDiff <- 10^-5
  NGen  <- 0
  MeanFreqOld <- rep(0, NGenPerBatch)
  
  # Loop over generations and calculate new allele frequencies
#  cat("Running simulation ... ")
  while(NGen <= MaxNGen & PDiff <= PDiffThresh){
    
    MeanFreqNew <- NULL
    for (i in 1:NGenPerBatch){
      SampleProbs <- AlleleFreq * PCoeff / (1 + AlleleFreq * (PCoeff - 1) )
      AlleleFreq  <- c(rbinom(length(SampleProbs), 2*N, SampleProbs) / (2*N), 
                       AlleleFreqStart)
      blnDrop     <- AlleleFreq == 0 | AlleleFreq == 1
      NDropOut    <- c(NDropOut, sum(blnDrop))
      AlleleFreq  <- AlleleFreq[!blnDrop]
      MeanFreqNew <- c(MeanFreqNew, mean(AlleleFreq))
      NrAlleles   <- c(NrAlleles, length(AlleleFreq))
    }
    MeanFreq <- c(MeanFreq, MeanFreqNew)
    PDiff    <- t.test(MeanFreqNew, MeanFreqOld)$p.value
    MeanFreqOld <- MeanFreqNew
  }
#  cat("done!\n")
  
  # Calculate probabilities for observed frequency counts
  AlleleFreqCount    <- table(AlleleFreq)
  AlleleFreqRelCount <- AlleleFreqCount / sum(AlleleFreqCount)
  AlleleFreqVals     <- as.numeric(names(AlleleFreqCount))
  FreqProbs <- sapply(FreqCountV, function(x) {
    AlleleFreqRelCount %*% dbinom(x, SampleSize, AlleleFreqVals)
  })
  FreqProbs <- c(FreqProbs, 
                 AlleleFreqRelCount %*% (1 - pbinom(max(FreqCountV) + 1, 
                                         SampleSize, AlleleFreqVals)))
  
  log(dmultinom(SFS, sum(SFS), FreqProbs))
}



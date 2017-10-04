##############################################
#
# General description:
#
#   The following function explores a grid of selection parameters. The 
#   parameters explored are the shape (a) and scale (sc) parameters of a gamma 
#   distribution of fitness values. The function simulates for each parameter
#   combination a distribution of allele frequencies, compares them to observed
#   frequencies and performs a regression a la Beaumont et al 2002 Genetics
#

# Input:
#
#     ObservedFreq = vector of observed allele frequencies
#     FitMeans     = vector of mean fitness values
#     FitVars      = vector of fitness variances
#     PopSize      = population size
#     NrSamples    = number of alleles sampled (see GenerateAlleleFreq)
#     NrRep        = number of alleles generated (see GenerateAlleleFreq)
#     NrGen        = number of generations (see GenerateAlleleFreq)
#     MaxFreq      = maximum allele frequency in observed data 
#                    (frequencies above this value are assumed to occur in 
#                     reference genome)
#     ProbV        = vector of probabilities (used for quantiles)
#     Epsilon      = maximum deviation for values to be included in regression
#     BreakV       = vector of frequency breaks (used for histogram)
#     SummaryType  = character ("none", "ks", "quant", or "hist") indicating
#                    what summary type is used to compare observed and 
#                    simulated frequencies

# Output:
#   
#     FitMeansRep  = replicated vector of mean fitness values 
#     FitVarsRep   = replicated vector of fitness variances
#     aVals        = vector of shape parameters used
#     scVals       = vector of scale parameters used
#     GridSummary  = vector or matrix that summarizes the simulated allele
#                    frequency for each parameter combination (1 column per 
#                    parameter combination)
#     ObsSummary   = vector summarizing the observed allele frequencies
#     DiffMean     = vector of mean differences between observed and simulated 
#                    summary  
#     DiffMax      = vector of max differences between observed and simulated 
#                    summary  
#     LMFit_var    = linear fit for fitness variance 
#     LMFit_mean   = linear fit for fitness means 
#     VarSample    = posterior sample of fitness variance 
#     MeanSample   = posterior sample of fitness mean 

# Comments:
#
#     Requires function GenerateAlleleFreq

##############################################

ExploreSelectionParameterGrid <- function(ObservedFreq, 
                                          FitMeans = seq(0.5, 1.5, 0.1),
                                          FitVars  = 0:10 + 0.1,
                                          PopSize  = 10^4,
                                          NrSamples = 10^4,
                                          NrRep = 10^4,
                                          NrGen = 10^3,
                                          MaxFreq = 0.5,
                                          ProbV = seq(0.05, 0.95, 0.1),
                                          Epsilon = 0.1,
                                          BreakV = seq(0, 1, 0.02),
                                          SummaryType = c("none", "ks", "quant", "hist")){
  
  # Repeat mean and variances of fitness values so that each gets combined with
  # the full range of variances
  FitMeansRep <- rep(FitMeans, length(FitVars))
  FitVarsRep  <- rep(FitVars, each = length(FitMeans))

  # Generate vectors of shape (aVals) and scale (scVals) parameter values
  aVals  <- FitMeansRep^2 / FitVarsRep
  scVals <- FitVarsRep / FitMeansRep
  
  # Replace frequencies above maxium value by 1
  ObservedFreq[ObservedFreq > MaxFreq] <- 1
  
  # Determine the type of summary
  SummaryType <- SummaryType[1]
  
  # Loop over grid values and calculate summaries
  GridSummary <- switch(SummaryType,
                        'none' = sapply(1:length(aVals), function(i){
                          GenerateAlleleFreq(aVals[i], scVals[i],  n = PopSize, 
                                             NrGen = NrGen, NrRep = NrRep,
                                             NrSamples = NrSamples)
                        }),
                        'ks'  = sapply(1:length(aVals), function(i){
                          # Simulate allele frequencies
                          AlleleFreq <- GenerateAlleleFreq(aVals[i], scVals[i],  n = PopSize, 
                                                           NrGen = NrGen, NrRep = NrRep,
                                                           NrSamples = NrSamples)
                          AlleleFreq[AlleleFreq>MaxFreq] <- 1
                          
                          # Calculate Kolmogorov-Smirnov statistic for the difference
                          cat("Calculating Kolmogorov-Smirnov test\n\n")
                          ks.test(AlleleFreq, ObservedFreq)$statistic
                        }),
                        'quant' = sapply(1:length(aVals), function(i){
                          # Simulate allele frequencies
                          AlleleFreq <- GenerateAlleleFreq(aVals[i], scVals[i],  n = PopSize, 
                                                           NrGen = NrGen, NrRep = NrRep,
                                                           NrSamples = NrSamples)
                          AlleleFreq[AlleleFreq > MaxFreq] <- 1
                          
                          # Calculate qunatiles for simulated frequencies
                          cat("Calculating quantiles\n\n")
                          quantile(AlleleFreq, ProbV)
                        }),
                        'hist' = sapply(1:length(aVals), function(i){
                          # Simulate allele frequencies
                          AlleleFreq <- GenerateAlleleFreq(aVals[i], scVals[i],  n = PopSize, 
                                                           NrGen = NrGen, NrRep = NrRep,
                                                           NrSamples = NrSamples)
                          AlleleFreq[AlleleFreq > MaxFreq] <- 1
                          
                          # Calculate qunatiles for simulated frequencies
                          cat("Calculating histogram\n\n")
                          Dens <- hist(AlleleFreq, BreakV, plot = F)$density
                          Dens / sum(Dens)
                        })
  )
  
  # Observed summary
  ObsSummary <- switch(SummaryType,
                       'none'  = 0,
                       'ks'    = 0,
                       'quant' = quantile(ObservedFreq, ProbV),
                       'hist'  = {
                         Dens <- hist(ObservedFreq, BreakV)$density
                         Dens / sum(Dens)}
  )
  
  # Estimate intercept via regression and generate a sample of means and
  # varuiances of the fitness distribution
  DiffMat    <- as.matrix(GridSummary) - ObsSummary
  AbsDiffMat <- abs(DiffMat) / mean(GridSummary)
  DiffMean   <- colMeans(AbsDiffMat)
  DiffMax    <- apply(AbsDiffMat, 2, max)
  blnE       <- DiffMax < Epsilon
  XMat       <- t(DiffMat)[blnE, ]
  if (sum(blnE) > 5 & SummaryType != 'none'){
    LMFit_var  <- lm(FitVarsRep[blnE]       ~ XMat)
    LMFit_mean <- lm(FitMeansRep[blnE] ~ XMat)
    VarSample <- XMat %*% LMFit_var$coefficients[-1]
    MeanSample <- XMat %*% LMFit_mean$coefficients[-1]
  } else {
    warning("Not enough values crossed threshold. No regression!\n")
    LMFit_var  <- NA
    LMFit_mean <- NA
    VarSample  <- NA
    MeanSample <- NA
    
  }
  
  # Return values in a list
  list(FitMeansRep  = FitMeansRep, FitVarsRep = FitVarsRep,
       aVals = aVals, scVals = scVals, 
       GridSummary  = GridSummary, ObsSummary = ObsSummary, DiffMean = DiffMean, 
       DiffMax = DiffMax,
       LMFit_var = LMFit_var, LMFit_mean = LMFit_mean, VarSample = VarSample,
       MeanSample = MeanSample)
}

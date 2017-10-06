##############################################
#
# General description:
#
#   The following function calculates the maximum difference between two 
#   cumulative distributions (as Kolmogorov-Smirnov but does not use
#   absolute values)
#

# Input:
#
#     Obs1 = vector of observed values
#     Obs2 = vector of observed values

# Output:
#   
#     Cum1  = cumulative distribution of obs1 
#     Cum2  = cumulative distribution of obs2 
#     idxMax = index of maximum difference between cumulative distributions
#     MaxDiff = maximum difference between cumulative distributions
## Comments:
#
#     Requires function GenerateAlleleFreq

##############################################

CumDistDiff <- function(Obs1, Obs2){
  PooledObs <- sort(c(Obs1, Obs2))
  Cum1 <- sapply(PooledObs, function(x) sum(Obs1 <= x)) / length(Obs1)
  Cum2 <- sapply(PooledObs, function(x) sum(Obs2 <= x)) / length(Obs2)
  CumDiff <- Cum1 - Cum2
  idxMax <- which.max(abs(CumDiff))
  list(PooledObs = PooledObs, Cum1 = Cum1, Cum2 = Cum2, idxMax = idxMax, 
       MaxDiff = CumDiff[idxMax])
}

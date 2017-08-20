##############################################
#
# General description:
#
#   The following function repeatedly samples a subset of observations,
#   calculates quantiles of the sampled subset and retains for each quantile
#   the lower and upper quantile amoong the samples

# Input:
#
#     Vals: vector of values to sample from
#     SubsetSize: size of sampled subset
#     NrSamples: number of subsets to be sampled
#     QuantV: vector of quantiles to pe sampled
#     LowerQ, UpperQ: lower and upper quantile of sample quantiles to be retained
#     

# Output:
#   
#    QMat: matrix of observed qunatiles (1st row), lower qunatile (2nd row) and
#          upper quantile (3rd row) for each sampled quantile
#    SampleMeans: vector of sampled means

##############################################

######                                      
# Source packages and set parameters  
######
SampleQuantiles <- function(Vals, SubsetSize, NrSamples = 1000, 
                            QuantV = seq(0, 1, 0.1), LowerQ = 0.025, UpperQ = 0.975){
  SampleMeans <- rep(NA, NrSamples)
  SampledQMat <- matrix(nrow = NrSamples, ncol = length(QuantV))
  for (i in 1:NrSamples){
    SampleVals <- sample(Vals, SubsetSize)
    SampleMeans[i]  <- mean(SampleVals)
    SampledQMat[i,] <- quantile(SampleVals, QuantV)
  }
  QMat <- matrix(nrow = 3, ncol = length(QuantV))
  QMat[1,] <- quantile(Vals, QuantV)
  QMat[2:3,] <- apply(SampledQMat, 2, 
                      FUN = function(x) quantile(x, c(LowerQ, UpperQ)))
  list(QMat = QMat, SampleMeans = SampleMeans)
}
  
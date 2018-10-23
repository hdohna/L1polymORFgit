# The following script explores functions used to estimate selection 
# coefficients

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Set parameters
N = 10^4
SampleSize = 2*2504
k = 5000
s = -0.001

# Calculate distributions of frequencies for different selection coefficients
SelectCoeffs <- seq(-0.001, 0.001, 0.0005)
Freqs <- seq(0.01, 0.99, 0.01)
DistMat1 <- sapply(SelectCoeffs, function(x){
  sapply(Freqs, function(y) {
    IntConst <- integrate(function(z) AlleleFreqTime(z, x, N), 0, 1)$value
    AlleleFreqTime(y, x, N) / IntConst
  })
})
Cols <- rainbow(ncol(DistMat1))
plot(Freqs, DistMat1[,1], col = Cols[1], type = "l")
for (i in 2:ncol(DistMat1)){
  lines(Freqs, DistMat1[,i], col = Cols[i])
}
legend("topright", legend = SelectCoeffs, lty = 1, col = Cols)

# Calculate count probabilities
CountVals <- seq(1, SampleSize, 10)
CountProbs1 <- sapply(SelectCoeffs, function(x){
  sapply(CountVals, function(y) {
    AlleleFreqSample(y, x, N, SampleSize = SampleSize, 
                     blnUseFPrime = F) 
  })
})
CountProbs2 <- sapply(SelectCoeffs, function(x){
  sapply(CountVals, function(y) {
    AlleleFreqSample(y, x, N, SampleSize = SampleSize, 
                     blnUseFPrime = T) 
  })
})
Cols <- rainbow(ncol(CountProbs1))
plot(CountVals, CountProbs1[,1], col = Cols[1], type = "l")
for (i in 2:ncol(CountProbs1)){
  lines(CountVals, CountProbs1[,i], col = Cols[i])
}
plot(CountVals, CountProbs2[,1], col = Cols[1], type = "l", xlim = c(1, 500))
for (i in 2:ncol(CountProbs2)){
  lines(CountVals, CountProbs2[,i], col = Cols[i])
}
plot(CountVals, exp(CountProbs2[,1]), col = Cols[1], type = "l")
for (i in 2:ncol(CountProbs2)){
  lines(CountVals, exp(CountProbs2[,i]), col = Cols[i])
}


LChooseVals <- sapply(CountVals, function(k) lchoose(SampleSize, k))
IntVals <- sapply(CountVals, function(k) {
  integrate(function(x) AlleleFreqTime(x, 0, N) * (dbinom(k, SampleSize, x) /
                                                     choose(SampleSize, k)), 0, 1)$value
  })
IntConstVal <- integrate(function(x) AlleleFreqTime(x, s, N), 
                         0, 1)$value - 
  integrate(function(x) (1 - x)^SampleSize * AlleleFreqTime(x, s, N), 
            0, 1)$value

plot(CountVals, LChooseVals, xlim = c(1, 300))
plot(CountVals, log(IntVals), xlim = c(1, 300))
plot(CountVals, LChooseVals + log(IntVals), xlim = c(1, 300))
sort(IntVals)

dbinom(500, SampleSize, 0.01)
ExactLL <- sapply(seq(1, 201, 10), function(x) AlleleFreqSample(x, -0.001, N, 
                                                                SampleSize = SampleSize, blnUseFPrime = T))

ApproxLLMat <- sapply(seq(0.01, 0.99, 0.01), function(y){
  sapply(seq(1, 201, 10), function(x) AlleleFreqSampleApprox(x, 
                                                             -0.001, N, SampleSize = SampleSize, blnUseFPrime = T, MaxW = y))
})
x = 1
sqDiff <- sapply(1:ncol(ApproxLLMat), function(x) mean((ExactLL[1:17] - ApproxLLMat[1:17,x])^2))
plot(sqDiff)
which.min(sqDiff)
plot(ExactLL, ApproxLLMat[,which.min(sqDiff)])

sVals <- seq(-0.002, 0.001, 0.0005)
CountVals <- c(seq(1, 171, 10), seq(4801, 5001, 10))
ExactLLMat <- sapply(sVals, function(y){
  sapply(CountVals, function(x) AlleleFreqSample(x, y, N, 
                                                 SampleSize = SampleSize, blnUseFPrime = T))
})
Cols <- rainbow(ncol(ExactLLMat))
plot(CountVals, ExactLLMat[,1], col = Cols[1], type = "l")
for (i in 2:ncol(ExactLLMat)){
  lines(CountVals, ExactLLMat[,i], col = Cols[i])
}


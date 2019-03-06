-AlleleFreqLogLik_4Par(
  Freqs = round(L1TotData$Freq[!blnNA], 0),
  Counts = rep(1, sum(!blnNA)),
  Predict = PredictMat[!blnNA, 1:3],
  a = -0.0001, b = 0, c = 0, d = 0, N = PopSize,
  SampleSize = L1TotData$SampleSize[!blnNA],
  blnIns = L1TotData$blnIns[!blnNA], 
  LogRegCoeff = LogRegL1Ref$coefficients,
  DetectProb = L1TotData$DetectProb[!blnNA])

Freqs = round(L1TotData$Freq[!blnNA], 0)
Counts = rep(1, sum(!blnNA))
Predict = PredictMat[!blnNA, 1:3]
a = 0
b = 0 
c = 0 
d = 0 
N = PopSize
SampleSize = L1TotData$SampleSize[!blnNA]
blnIns = L1TotData$blnIns[!blnNA]
LogRegCoeff = LogRegL1Ref$coefficients
DetectProb = L1TotData$DetectProb[!blnNA]

sVals <- as.matrix(Predict) %*% c(b, c, d)
i = 1
Counts[i] * AlleleFreqSample(Freqs[i], 
                             a + sVals[i], N, 
                             SampleSize = SampleSize[i], 
                             blnIns = blnIns[i],
                             LogRegCoeff = LogRegCoeff, 
                             DetectProb = DetectProb[i])


k = Freqs[i]
s = a + sVals[i]
N = 10^5 
SampleSize = SampleSize[i] 
blnIns = blnIns[i]
LogRegCoeff = LogRegCoeff
DetectProb = DetectProb[i]
blnIns = blnIns[i]

ProbRef <- function(x) {
  ExpCoeff <- exp(LogRegCoeff[1] + LogRegCoeff[2] * x)
  ExpCoeff / (1 + ExpCoeff)
}
plot()
ProbRef(seq(0, 1, 0.01))

AlleleFreqTime(10^-3, s, N)
AlleleFreqTime(10^-5, s, N)
AlleleFreqTime(10^-6, s, N)
AlleleFreqTime(10^-10, s, N)

SampleSize = 2000
OffSet <- 10^-4
IntConst <- integrate(function(x) {
  AlleleFreqTime(x, s, N) * (1 - ProbRef(x))
},  OffSet, 1 - OffSet)$value - 
  integrate(function(x) {
    AlleleFreqTime(x, s, N) * (1 - DetectProb * x)^SampleSize * (1 - ProbRef(x))
  },  OffSet, 1 - OffSet)$value


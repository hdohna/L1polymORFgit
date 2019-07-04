library(expint)


# The equation below is from Bossinoit et al 2006 for online integration tools

# x < 1/2N
(e^(-s)*(e^s - e^(2*N*s))*(e^(2*N*s*x) - 1))/((e^(2*N*s) - 1)*N*s*(x - 1)*x) 

# Result of integral from 0 to 1/(2n) (https://www.integral-calculator.com/)
(e^(-s)*(e^(2*n*s)-e^s)*(ln(abs(x-1)/abs(x))+expintegral_ei(2*n*s*x)-e^(2*n*s)*expintegral_ei(2*n*s*(x-1))))/(n*s*(e^(2*n*s)-1))

# Result of integral from 0 to 1/(2n) multiplied with ((https://www.integral-calculator.com/)
(e^(-s)*(e^(2*n*s)-e^s)*(ln(abs(x-1)/abs(x))+expintegral_ei(2*n*s*x)-e^(2*n*s)*expintegral_ei(2*n*s*(x-1))))/(n*s*(e^(2*n*s)-1))

# First derivative of numerator with respect to s:
# https://www.derivative-calculator.net/
-e^(-s)*(((2*n*x + 2*n - 1)*e^(2*n*s) - 2*n*x*e^s) * e^(2*n*x*s) + (1 - 2*n)*e^(2*n*s))
-1*((2*n - 1) + (1 - 2*n)) # s = 0
0 # s = 0 (verified once)

# Second derivative of numerator
-e^(-s)*(((4*n^2*x^2+(8*n^2-4*n)*x+4*n^2-4*n+1)*e^(2*n*s)-4*n^2*x^2*e^s)*e^(2*n*x*s)+(-4*n^2+4*n-1)*e^(2*n*s))
-1*(((4*n^2*x^2 + (8*n^2 - 4*n)*x + 4*n^2 - 4*n + 1) - 4*n^2*x^2) + (-4*n^2 + 4*n - 1)) # s = 0
-1*((8*n^2 - 4*n)*x) # s = 0
4*x*n*(1 - 2*n) # s= 0 (verified once)

# First derivative of numerator with respect to s:
n*(x-1)*x*((2*n*s+1)*e^(2*n*s)-1) 
n*(x-1)*x*(1-1) = 0 # s=0

# Second derivative of denominator
4*n^2*(x-1)*x*(n*s+1)*e^(2*n*s)
4*n^2*(x-1)*x# s = 0

# Ratio of second derivatives for s = 0
4*x*n*(1 - 2*n) / (4*n^2*x*(x - 1)) = (1 - 2*n) / (n*(x - 1))
    
    
# x > 1/2N
(-e^(-s)*(e^s - 1)*(e^(2*N*s) - e^(2*N*s*x)))/((e^(2*N*s) - 1)*N*s*(x - 1)*x) 

# Result (https://www.integral-calculator.com/)
(e^(-s)*(e^s-1)*(e^(2*n*s)*(ln(abs(x))+expintegral_ei(2*n*s*(x-1))-ln(abs(x-1)))-expintegral_ei(2*n*s*x)))/(n*s*(e^(n*s)-1)*(e^(n*s)+1))

# First derivative of numerator with respect to s:
e^(-s)*((2*n*x*e^s-2*n*x+1)*e^(2*n*x*s)+(-2*n*e^s+2*n-1)*e^(2*n*s))
1*((2*n*x - 2*n*x + 1) + (-2*n + 2*n - 1)) # s = 0
1*(1 - 1) = 0 # s=0

# Second derivative of numerator with respect to s
e^(-s)*((4*n^2*x^2*e^s-4*n^2*x^2+4*n*x-1)*e^(2*n*x*s)+(-4*n^2*e^s+4*n^2-4*n+1)*e^(2*n*s))
(4*n^2*x^2 - 4*n^2*x^2 + 4*n*x - 1) + (-4*n^2 + 4*n^2 - 4*n + 1) # s = 0
4*n*(x - 1) # s = 0 (verified once)

# Ratio of second derivatives for s = 0
4*n*(x - 1) / (4*n^2*x*(x - 1)) = 1 / (n*x)

# Checking that ratios are the same when x = 1/(2n)
(1 - 2*n) / (n*(x - 1)) = 1 / (n*x) # plugging in x = 1/(2n)
(1 - 2*n) / (1/2 - n)   = 1 / (1/2)
2*(1/2 - n) / (1/2 - n) = 2 = 1 / (1/2)


# R function for distribution of frequency as function of selection coefficient
FreqD <- function(y, s, N) {
  exp(-s)*(((1 + sign(1/(2*N) - y)) / 2 * (exp(s) - exp(2*N*s))*(exp(2*N*s*y) - 1))
   - (1 + sign(y - 1/(2*N))) / 2  * ((exp(s) - 1)*(exp(2*N*s) - exp(2*N*s*y)))) / 
     ((exp(2*N*s) - 1)*N*s*(y - 1)*y)  /
    ((exp(-s)*(exp(2*N*s)-exp(s))*(log(abs(y-1)/abs(y))+expint_Ei(2*N*s*y)-exp(2*N*s)*expint_Ei(2*N*s*(y-1))))/(N*s*(exp(2*N*s)-1))
 + (exp(-s)*(exp(s)-1)*(exp(2*N*s)*(log(abs(y))+expint_Ei(2*N*s*(y-1))-log(abs(y-1)))-expint_Ei(2*N*s*y)))/(N*s*(exp(N*s)-1)*(exp(N*s)+1)))
  
}
FreqD(y = 0.1, s = 0.1, N = 10^4)
yVals <- seq(0, 1, 0.01)
Fvals1 <- sapply(yVals, function(x) FreqD(x, s = 0.02, N = 10^4))
Fvals2 <- sapply(yVals, function(x) FreqD(x, s = 0.001, N = 10^4))
Fvals3 <- sapply(yVals, function(x) FreqD(x, s = 10^(-4), N = 10^4))
Fvals4 <- sapply(yVals, function(x) FreqD(x, s = 10^(-5), N = 10^4))
Fvals5 <- sapply(yVals, function(x) FreqD(x, s = -0.001, N = 10^4))
plot(yVals, Fvals1, type = "l", col = "blue", xlim = c(0, 0.2))
lines(yVals, Fvals2, col = "green")
lines(yVals, Fvals3, col = "red")
lines(yVals, Fvals4, col = "orange")
lines(yVals, Fvals5, col = "brown")

FreqD_unst <- function(y, s, N) {
  exp(-s)*(((1 + sign(1/(2*N) - y)) / 2 * (exp(s) - exp(2*N*s))*(exp(2*N*s*y) - 1))
           - (1 + sign(y - 1/(2*N))) / 2  * ((exp(s) - 1)*(exp(2*N*s) - exp(2*N*s*y)))) / 
    ((exp(2*N*s) - 1)*N*s*(y - 1)*y)
  }
FreqD_unst(y = 0.1, s = 0.05, N = 10^4)
MeanY <- function(yVals, s, N){
  FreqVals <- sapply(yVals, function(y) y * FreqD_unst(y, s, N))
  if (any(is.na(FreqVals))){
    warning(sum(is.na(FreqVals)), " frequency values are NA\na")
  }
  sum(FreqVals[!is.na(FreqVals)]) / 
    sum(yVals[!is.na(FreqVals)])
}

# Plot mean frequency per selection value
SelVals <- seq(-0.03, 0.03, 0.001)
YperSel <- sapply(SelVals, function(x){
  MeanY(y = seq(0, 1, 0.01), s = x, N = 10^4)
})
plot(SelVals, YperSel, type = "l") 
plot(SelVals, log10(YperSel), type = "l") 

FreqD_unst <- function(y, s, N) {
  exp(-s)*(((1 + sign(1/(2*N) - y)) / 2 * (exp(s) - exp(2*N*s))*(exp(2*N*s*y) - 1))
           - (1 + sign(y - 1/(2*N))) / 2  * ((exp(s) - 1)*(exp(2*N*s) - exp(2*N*s*y)))) / 
    ((exp(2*N*s) - 1)*N*s*(y - 1)*y)
}

FreqD_st <- function(y, s, N){
  FreqD_unst(y, s, N) / integrate(function(x) FreqD_unst(x, s, N), 0, 1)$value
}
FreqD_sample2 <- function(y, s, N, SampleSize = 10^3){
  (1 - (1 - y)^SampleSize) * FreqD_st(y, s, N) / 
    integrate(function(x) (1 - (1 - x)^SampleSize) * FreqD_st(x, s, N), 0, 1)$value
}
FreqD_sample <- function(y, s, N, SampleSize = 10^3){
  (1 - (1 - y)^SampleSize) * FreqD_unst(y, s, N) / 
    integrate(function(x) (1 - (1 - x)^SampleSize) * FreqD_unst(x, s, N), 0, 1)$value
}
FreqD_st(0, 0.01, 10^4)
FreqD_unst(0, 0.01, 10^4)
Fvals_st1 <- sapply(yVals, function(x) FreqD_sample(x, s = 0.02, N = 10^4))
Fvals_st2 <- sapply(yVals, function(x) FreqD_sample(x, s = 0.001, N = 10^4))
Fvals_st3 <- sapply(yVals, function(x) FreqD_sample(x, s = 10^(-4), N = 10^4))
Fvals_st4 <- sapply(yVals, function(x) FreqD_sample(x, s = 10^(-5), N = 10^4))
Fvals_st5 <- sapply(yVals, function(x) FreqD_sample(x, s = 0, N = 10^4))
Fvals_st5 <- sapply(yVals, function(x) FreqD_sample(x, s = -0.001, N = 10^4))
plot(yVals, Fvals_st1, type = "l", col = "blue")
lines(yVals, Fvals_st2, col = "green")
lines(yVals, Fvals_st3, col = "red")
lines(yVals, Fvals_st4, col = "orange")
lines(yVals, Fvals_st5, col = "brown")

FreqD_sample2(0.1, s = 0.001, N = 10^4)
FreqD_sample(0.1, s = 0.001, N = 10^4)

# Log-likelihood 
LogLik <- function(Freqs, Counts, s, N){
  if (length(Freqs) != length(Counts)){
    stop("Freqs and Counts vector have to have the same length\n")
  }
  SampleSize <- sum(Counts)
  sum(sapply(1:length(Freqs), function(i){
    Counts[i] * log(FreqD_sample(Freqs[i]/SampleSize, s, N, 
                                 SampleSize = SampleSize))
  }))
  
}


LogLik(Freqs = c(1, 2), Counts = c(100, 10), s = 10^(-4), N = 10^4)

FreqD_sample(1/110, s = 10^(-4), N = 10^4, 
             SampleSize = 110)

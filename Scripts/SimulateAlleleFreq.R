# The following script simulates allele frequencies for different selection
# coefficients

##########################################
#                                        #
#           Set parameters               #
#                                        #
##########################################

# Set population size
n = 10^4

# Set selection coefficients
SelCoeff <- seq(-0.1, 0.1, 0.01)

# Set the number of replicates per selection coefficient
NrRep <- 1000

# Set the number of generations
NrGen <- 10000

##########################################
#                                        #
#           Simulate                     #
#                                        #
##########################################

# Create a vector of replicated selection coefficients
SelCoeffRep <- rep(SelCoeff, NrRep)
PCoeff      <- 1 + SelCoeffRep
cat("Length of allele frequency vector is", length(SelCoeffRep), "\n")

# Initialize vector of allele frequencies
AlleleFreq <- rep(1/(2*n), length(SelCoeffRep))

# Loop over generations and calculate new allele frequencies
cat("Running simulation ... ")
for (i in 1:NrGen){
  SampleProbs <- AlleleFreq * PCoeff / (1 + AlleleFreq * (PCoeff - 1) )
  AlleleFreq  <- rbinom(length(SampleProbs), 2*n, SampleProbs) / (2*n)
  AlleleFreq  <- pmin(1 - 1/(2*n), AlleleFreq)
  AlleleFreq  <- pmax(1/(2*n), AlleleFreq)
}
cat("done!\n")
hist(AlleleFreq, seq(0, 1, 1/n))
hist(AlleleFreq[SelCoeff == 0], seq(0, 1, 1/n))
hist(AlleleFreq[abs(SelCoeff - 0.01) < 10^-10], seq(0, 1, 1/n))
hist(AlleleFreq[abs(SelCoeff + 0.01) < 10^-10], seq(0, 1, 1/n))

##########################################
#                                        #
#           Alternative simulation                  #
#                                        #
##########################################

# Number to be added per generation
NAdd <- 100

# Initialize vector of allele frequencies
PCoeff <- 1 + 10^-3
AlleleFreqStart <- rep(1/(2*n), NAdd)
AlleleFreq      <- AlleleFreqStart
NDropOut        <- 0
MeanFreq        <- 0
NrAlleles       <- length(AlleleFreq)
# Loop over generations and calculate new allele frequencies
cat("Running simulation ... ")
for (i in 1:10000){
  SampleProbs <- AlleleFreq * PCoeff / (1 + AlleleFreq * (PCoeff - 1) )
  AlleleFreq  <- c(rbinom(length(SampleProbs), 2*n, SampleProbs) / (2*n), AlleleFreqStart)
  blnDrop     <- AlleleFreq == 0 | AlleleFreq == 1
  NDropOut    <- c(NDropOut, sum(blnDrop))
  AlleleFreq  <- AlleleFreq[!blnDrop]
  MeanFreq    <- c(MeanFreq, mean(AlleleFreq))
  NrAlleles    <- c(NrAlleles, length(AlleleFreq))
}
cat("done!\n")
# hist(AlleleFreq, seq(0, 1, 1/n))
# hist(AlleleFreq[SelCoeff == 0], seq(0, 1, 1/n))
# hist(AlleleFreq[abs(SelCoeff - 0.01) < 10^-10], seq(0, 1, 1/n))
# hist(AlleleFreq[abs(SelCoeff + 0.01) < 10^-10], seq(0, 1, 1/n))
# plot(NDropOut, type = "l")
# plot(MeanFreq, type = "l")
plot(NrAlleles, type = "l")
# dbinom(1, 5000, AlleleFreq)
plot(sapply(1:50, function(x) sum(dbinom(x, 5000, AlleleFreq))))

# Function to calculate the mean frequency above a particular threshold 
# frequency, given the mean and variance of a normal selection
# coefficient distribution
MeanFreqNorm <- function(Mu = 0, SD = 0.01, Thresh = 0.5){
  Distn   <- dnorm(SelCoeffRep, Mu, sd = SD)
  Distn[Distn == Inf] <- 10^6
  RelFreq <- Distn * AlleleFreq / sum(Distn)
  mean(RelFreq[AlleleFreq >= Thresh])
}
Distn   <- dnorm(SelCoeffRep, -0.01, 0)
sum(Distn == Inf)
sum(SelCoeffRep == -0.01)

# Plot the ratio of mean frequency above 0.5 for high vs low sd of selection
# coefficient
cat("Calculating frequency ratios for different selection variances ... ")
StDevV   <- seq(10^(-6), 0.1, 0.001)
NegSelCoeff <- -1.9/n
Threshhold <- 0.9
MeanF0   <- sapply(StDevV, function(x) MeanFreqNorm(SD = x, Thresh = Threshhold))
MeanF1   <- sapply(StDevV, function(x) MeanFreqNorm(Mu = NegSelCoeff, SD = x, 
                                                    Thresh = Threshhold))
MaxRatio <- max(MeanF1 / MeanF1[2], na.rm = T)
plot(StDevV,  MeanF0 / MeanF0[2], type = "l", ylim = c(0, MaxRatio),
     xlab = "Standard deviation", ylab = "Ratio in mean frequency")
lines(StDevV, MeanF1 / MeanF1[2], lty = 2)
cat("done!\n")

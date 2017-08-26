# The following script fits parameters of a distribution of fitness values
# to a distribution of allele frequencies

##########################################
#                                        #
#           Set parameters               #
#                                        #
##########################################

# Set population size
n = 10^4

# Set the number of replicates 
NrRep <- 10000

# Set the number of generations
NrGen <- 1000

# Parameters of the gamma distribution of selection coefficients
Gshape = 90
GSscale = 1/100

# Plot gamma distn for comparison
xVals <- seq(0, 2, 0.01)
plot(xVals, dgamma(xVals, shape = Gshape, scale = GSscale))


##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################
cat("Loading and processing data ... ")

# Read in data on transposable elements
MRIP <- read.delim("D:/L1polymORF/Data/eul1db_MRIP.txt", skip = 5)
SRIP <- read.delim("D:/L1polymORF/Data/eul1db_SRIP.txt", skip = 5)
IDmatch        <- match(MRIP$X.mrip_accession.no, SRIP$mrip)
MRIP$integrity <- SRIP$integrity[IDmatch]
MRIP$subgroup  <- SRIP$sub_group[IDmatch]
MRIP <- MRIP[MRIP$subgroup == "L1-Ta", ]
cat("done!\n")

##########################################
#                                        #
#           Simulate                     #
#                                        #
##########################################

# Function to generate allele frequencies based on parameters of a distribution
# of Selection coefficients
GenerateAlleleFreq <- function(Gshape, GSscale, n = 10^4, NrGen = 10^3,
                               NrRep = 10^4, NrSamples = 10^3){
  cat("Simulating allele frequencies with Gshape =", Gshape, 
      "and GSscale = ", GSscale, " ....")
  # Create a vector of replicated selection coefficients
  PCoeff  <- rgamma(NrRep, shape = Gshape, scale = GSscale)
  
  # Initialize vector of allele frequencies
  AlleleFreq <- rep(1/(2*n), length(PCoeff))
  
  # Loop over generations and calculate new allele frequencies
  for (i in 1:NrGen){
    if (any(is.na(AlleleFreq))) browser()
    SampleProbs <- AlleleFreq * PCoeff / (1 + AlleleFreq * (PCoeff - 1) )
    AlleleFreq  <- rbinom(length(SampleProbs), 2*n, SampleProbs) / (2*n)
    AlleleFreq  <- pmin(1 - 1/(2*n), AlleleFreq)
    AlleleFreq  <- pmax(1/(2*n), AlleleFreq)
  }
  cat("done!\n")
  sample(AlleleFreq, NrSamples, prob = AlleleFreq)
}

# Function to calculate Kolmogorov-Smirnov Statistic for the difference between
# simulated and observed frequencies.
DiffAlleleFreqKS <- function(ObservedFreq, Gshape, GSscale, n = 10^4, NrGen = 10^3,
                             NrRep = 10^4, NrSamples = 10^3){ 
   
  # Simulate allele frequencies
  AlleleFreq <- GenerateAlleleFreq(Gshape, GSscale, n = n, NrGen = NrGen,
                                 NrRep = NrRep, NrSamples = NrSamples)
  
  # Calculate Kolmogorov-Smirnov statistic for the difference
  cat("Calculating Kolmogorov-Smirnov test\n\n")
  ks.test(AlleleFreq, ObservedFreq)$statistic
}
DiffAlleleFreqKS(MRIP$pseudoallelefreq,70, 1/100, 100)
DiffAlleleFreqKS(MRIP$pseudoallelefreq,10^(3), 10^(-3), 10)

# Find optimal shape, scale and population size
OptParN1000 <- optim(par = c(70, 1/100), fn = function(x) {
  DiffAlleleFreqKS(MRIP$pseudoallelefreq, Gshape = x[1], GSscale = x[2], n = 1000)},
  lower = c(10^(-3), 10^(-9), 10)
)
OptParN100 <- optim(par = c(70, 1/100), fn = function(x) {
  DiffAlleleFreqKS(MRIP$pseudoallelefreq, Gshape = x[1], GSscale = x[2], n = 100)},
  lower = c(10^(-3), 10^(-6), 10)
)

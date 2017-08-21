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
GenerateAlleleFreq <- function(Gshape, GSscale, n = 10^4, NrGen = 10^3){
  
  # Create a vector of replicated selection coefficients
  PCoeff  <- rgamma(NrRep, shape = Gshape, scale = GSscale)
  
  # Initialize vector of allele frequencies
  AlleleFreq <- rep(1/(2*n), length(PCoeff))
  
  # Loop over generations and calculate new allele frequencies
  for (i in 1:NrGen){
    SampleProbs <- AlleleFreq * PCoeff / (1 + AlleleFreq * (PCoeff - 1) )
    AlleleFreq  <- rbinom(length(SampleProbs), 2*n, SampleProbs) / (2*n)
    AlleleFreq  <- pmin(1 - 1/(2*n), AlleleFreq)
    AlleleFreq  <- pmax(1/(2*n), AlleleFreq)
  }
  AlleleFreq
}

AlleleFreq <- GenerateAlleleFreq(10000, 1/10000, NrGen = 1000)
hist(AlleleFreq, breaks = seq(0, 1, 0.001))
hist(MRIP$pseudoallelefreq)
ks.test(AlleleFreq, MRIP$pseudoallelefreq)$statistic

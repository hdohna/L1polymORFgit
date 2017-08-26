# The following script fits parameters of a distribution of fitness values
# to a distribution of allele frequencies using analytical methods

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
#     Calculate likelihood               #
#                                        #
##########################################

# Set parameters
a <- 7
b  <- 10
y  <- 0.01

# -(exp(-x)*(e^x-1)*(e^(2*n*x)-e^(2*n*y*x)))/(n*(y-1)*y*x*(e^(2*n*x)-1))
fa <- function(x) {
  b^a / gamma(a) * x^(a -1) * exp(-b * (x - 1)) * 
  (-exp(-(x - 1))*(exp(x - 1) - 1) * (exp(2*n*(x - 1)) - exp(2*n*y*(x - 1))) /
     (n*(y-1)*y*(x - 1)*(exp(2*n*(x - 1)) - 1)))
}
x <- 1.04
b^a / gamma(a) * x^(a -1) * exp(-b * (x - 1))* (-exp(-(x - 1))*(exp(x - 1) - 1)*(exp(2*n*(x - 1)) - exp(2*n*y*(x - 1))) /(n*(y-1)*y*(x - 1)*(exp(2*n*(x - 1))-1)))
b^a / gamma(a) * x^(a -1) * exp(-b * (x - 1))* (exp(-(x - 1))*(exp((x - 1))-1)*(exp(2*n*(x - 1))-exp(2*n*y*(x - 1))))/(n*(y-1)*y*(x - 1)*(exp(2*n*(x - 1))-1))
sVals <- c(1 - seq(10^(-10), 10^(-2), 10^(-5)), 1 + seq(10^(-10), 10^(-2), 10^(-5)))
plot(sVals, fa(sVals))
fa(1.0001)
integrate(function(x) fa(x), lower = 0, upper = Inf)

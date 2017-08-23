# The following script fits parameters of a distribution of fitness values
# to a distribution of allele frequencies, using the methods described in 
# Beaumont et al. 2002 Genetics

##########################################
#                                        #
#           Load packages                #
#                                        #
##########################################

# library(locfit)

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
MRIP <- read.delim("/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/eul1db_MRIP.txt", skip = 5)
SRIP <- read.delim("/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/eul1db_SRIP.txt", skip = 5)
IDmatch        <- match(MRIP$X.mrip_accession.no, SRIP$mrip)
MRIP$integrity <- SRIP$integrity[IDmatch]
MRIP$subgroup  <- SRIP$sub_group[IDmatch]
MRIP <- MRIP[MRIP$subgroup == "L1-Ta", ]
cat("done!\n")

# Test whether full-length and fragment L1 have different frequency
idxFull  <- which(MRIP$integrity == "full-length")
blnFragm <- MRIP$integrity %in% c("5prime-truncated", "3prime-truncated")
ks.test(MRIP$pseudoallelefreq[idxFull], MRIP$pseudoallelefreq[blnFragm])
var(MRIP$pseudoallelefreq[idxFull])
var(MRIP$pseudoallelefreq[blnFragm])
mean(MRIP$pseudoallelefreq[idxFull])
mean(MRIP$pseudoallelefreq[blnFragm])
hist(MRIP$pseudoallelefreq[idxFull], breaks = seq(0, 0.1, 0.0005))
hist(MRIP$pseudoallelefreq[blnFragm], breaks = seq(0, 0.1, 0.0005))

###################################################
#                                                 #
#  Functions to simulate allele frequencies       #
#                                                 #
###################################################

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

# Function to explore a grid of alpha values and fitness values
# Grid of parameter values
ExploreGrid <- function(ObservedFreq, 
                        aValsBasic = seq(1, 101, 10),
                        proPs      = seq(0.5, 1.5, 0.1),
                        PopSize  = 10^4,
                        blnPlot = T){
  
  # Repeat proportion (fitness) values so that each alpha gets combined with
  # the full range of fitness values
  proPsRep   <- rep(proPs, length(aValsBasic))
  aVals      <- rep(aValsBasic, each = length(proPs))
  bVals      <- 1/aVals * proPsRep
  
  # Evaluate difference according to Kolmogorov-Smirnov statistic 
  DiffKS <- sapply(1:length(aVals), function(i){
    DiffAlleleFreqKS(ObservedFreq, aVals[i], bVals[i], n = PopSize)
  })
  
  if (blnPlot){
    par(mfrow = c(2, 2))
    # Plot difference vs alpha
    Cols <- rainbow(length(proPs))
    plot(aVals, DiffKS, xlab = "alpha")
    for (i in 1:length(Cols)){
      blnProps <- proPsRep == proPs[i]
      points(aVals[blnProps], DiffKS[blnProps], col = Cols[i])
    }
    
    # Plot difference vs selection coefficient
    Cols <- rainbow(length(aValsBasic))
    plot(proPsRep, DiffKS, xlab = "Mean fitness")
    for (i in 1:length(Cols)){
      blnA <- aVals == aValsBasic[i]
      points(proPsRep[blnA], DiffKS[blnA], col = Cols[i])
    }
    
  }
  
  
  # Return values in a list
  list(proPsRep  = proPsRep, aVals = aVals, bVals = bVals, DiffKS = DiffKS)
}

###################################################
#                                                 #
#  Explore fit for simulations                    #
#                                                 #
###################################################

#######
# Analyze coarse grid
#######

cat("******  Exploring coarse grid    ***********\n\n")
# Results for distribution of full-length L1
ResultList1Full <- ExploreGrid(MRIP$pseudoallelefreq[idxFull],
                           aValsBasic = seq(1, 101, 10),
                           proPs = seq(0.2, 1.5, 0.1))

# Results for distribution of fragment L1
ResultList1Fragm <- ExploreGrid(MRIP$pseudoallelefreq[blnFragm],
                               aValsBasic = seq(seq(1, 101, 10)),
                               proPs = seq(0.2, 1.5, 0.1))

#######
# Analyze finer grid
#######

cat("******  Exploring fine grid    ***********\n\n")
# Results for distribution of full-length L1
ResultList2Full <- ExploreGrid(MRIP$pseudoallelefreq[idxFull],
                           aValsBasic = seq(70, 100, 2),
                   proPs = seq(0.8, 0.9, 0.001))
dev.copy2pdf("/srv/gsfs0/projects/levinson/hzudohna/L1InsertionLocation/L1FullFitDistnPlot.pdf")

# Results for distribution of fragment L1
ResultList2Full <- ExploreGrid(MRIP$pseudoallelefreq[blnFragm],
                               aValsBasic = seq(70, 100, 2),
                               proPs = seq(0.7, 0.9, 0.001))
dev.copy2pdf("/srv/gsfs0/projects/levinson/hzudohna/L1InsertionLocation/L1FragmFitDistnPlot.pdf")

#########
# Save results
#########

save.image("/srv/gsfs0/projects/levinson/hzudohna/L1InsertionLocation/SelectionParameterFit.RData")
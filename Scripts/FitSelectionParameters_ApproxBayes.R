# The following script fits parameters of a distribution of fitness values
# to a distribution of allele frequencies, using the methods described in 
# Beaumont et al. 2002 Genetics

##########################################
#                                        #
#           Load packages                #
#                                        #
##########################################

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# library(locfit)

##########################################
#                                        #
#           Set parameters               #
#                                        #
##########################################

# Path to ouput file
OutPath <- "D:/L1polymORF/Data/SelectionParameterFit_1000G_maxF0.3_FineGrid.RData"

# Set parameters for simulating allele frequencies
PopSize   <- 10^4
NrSamples <- 5000
NrRep     <- 10000
NrGen     <- 1000
MaxF      <- 0.3
PV        <- seq(0, 0.40, 0.1)
Epsilon   <- 0.01
BV        <- seq(0, 1, 0.02)

# mean and variance of fitnes values for grid
FitMeanGrid_Full <- seq(0.94, 0.98, 0.001)
FitVarGrid_Full  <- seq(0.001, 0.2, 0.01)
FitMeanGrid_Fragm <- seq(0.75, 0.85, 0.005)
FitVarGrid_Fragm  <- seq(1, 5, 0.1)

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################
cat("Loading and processing data ... ")

# Load L1 catalog GenomicRanges
load("D:/L1polymORF/Data/L1CatalogGRanges.RData")
load('D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData')

# Read in info about reference L1
L1Ref <- read.csv("D:/L1polymORF/Data/repeatsHg19_L1HS.csv")
L1Ref$Width  <- L1Ref$genoEnd - L1Ref$genoStart
L1RefNrFull  <- sum(L1Ref$Width >  6000, na.rm = T)
L1RefNrFragm <- sum(L1Ref$Width <= 5900, na.rm = T)

# Get frequency of full-length and fragment L1 in 1000 genome data
FreqFull_1000G  <- L1_1000G_reduced$Frequency[which(L1_1000G_reduced$InsLength > 6000)]
FreqFragm_1000G <- L1_1000G_reduced$Frequency[which(L1_1000G_reduced$InsLength <= 5900)]

# Add frequency of one for all insertions above maximum frequency and reference 
# insertions
NrAboveFull     <- sum(FreqFull_1000G >= MaxF)  + L1RefNrFull
NrAboveFragm    <- sum(FreqFragm_1000G >= MaxF) + L1RefNrFragm
FreqFull_1000G  <- c(FreqFull_1000G[FreqFull_1000G < MaxF], rep(1, NrAboveFull))
FreqFragm_1000G <- c(FreqFragm_1000G[FreqFragm_1000G < MaxF], 
                     rep(1, NrAboveFragm))

cat("done!\n")


###################################################
#                                                 #
#  Explore fit for simulations                    #
#                                                 #
###################################################


cat("******  Exploring  grid    ***********\n\n")
# Results for distribution of full-length L1
ResultListFull_1000G_hist <- ExploreSelectionParameterGrid(FreqFull_1000G, 
                              FitMeans = FitMeanGrid_Full,
                              FitVars  = FitVarGrid_Full,
                              PopSize  = PopSize,
                              NrSamples = NrSamples,
                              NrRep = NrRep,
                              NrGen = NrGen,
                              MaxFreq = MaxF,
                              ProbV = PV,
                              Epsilon = Epsilon,
                              BreakV = BV,
                              SummaryType = "hist")

# Results for distribution of full-length L1
# ResultListFull_1000G_quant <- ExploreSelectionParameterGrid(FreqFull_1000G, 
#                                 FitMeans = FitMeanGrid_Full,
#                                 FitVars  = FitVarGrid_Full,
#                                 PopSize  = PopSize,
#                                 NrSamples = NrSamples,
#                                 NrRep = NrRep,
#                                 NrGen = NrGen,
#                                 MaxFreq = MaxF,
#                                 ProbV = PV,
#                                 Epsilon = Epsilon,
#                                 BreakV = BV,
#                                 SummaryType = "quant")

# Results for distribution of fragment L1
ResultListFragm_1000G_hist <- ExploreSelectionParameterGrid(FreqFragm_1000G, 
                               FitMeans = FitMeanGrid_Fragm,
                               FitVars  = FitVarGrid_Fragm,
                               PopSize  = PopSize,
                               NrSamples = NrSamples,
                               NrRep = NrRep,
                               NrGen = NrGen,
                               MaxFreq = MaxF,
                               ProbV = PV,
                               Epsilon = Epsilon,
                               BreakV = BV,
                               SummaryType = "hist")

# Results for distribution of full-length L1
# ResultListFragm_1000G_quant <- ExploreSelectionParameterGrid(FreqFragm_1000G, 
#                                 FitMeans = FitMeanGrid_Fragm,
#                                 FitVars  = FitVarGrid_Fragm,
#                                 PopSize  = PopSize,
#                                 NrSamples = NrSamples,
#                                 NrRep = NrRep,
#                                 NrGen = NrGen,
#                                 MaxFreq = MaxF,
#                                 ProbV = PV,
#                                 Epsilon = Epsilon,
#                                 BreakV = BV,
#                                 SummaryType = "quant")
#########
# Save results
#########

save.image(OutPath)

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
OutPath <- "D:/L1polymORF/Data/SelectionParameterFit_Bossinot_KS_FineGrid.RData"

# Set parameters for simulating allele frequencies
PopSize   <- 10^4
NrSamples <- 5000
NrRep     <- 10000
NrGen     <- 1000
MaxF      <- 1
PV        <- seq(0, 0.40, 0.1)
Epsilon   <- 0.01
BV        <- seq(0, 1, 0.02)

# mean and variance of fitnes values for grid
FitMeanGrid_Full <- seq(0.9, 1.1, 0.01)
FitVarGrid_Full  <- seq(0.01, 0.51, 0.01)
FitMeanGrid_Fragm <- seq(0.9, 1.1, 0.01)
FitVarGrid_Fragm  <- seq(0.01, 0.51, 0.01)

# Specify summary type
SummaryType <- "ks"

# DataSource (1000G or Bossinot2006)
DataSource <- "Bossinot2006"

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################
cat("Loading and processing data ... ")

if (DataSource == "1000G"){
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
  FreqFull  <- c(FreqFull_1000G[FreqFull_1000G < MaxF], rep(1, NrAboveFull))
  FreqFragm <- c(FreqFragm_1000G[FreqFragm_1000G < MaxF], rep(1, NrAboveFragm))
  
} else {
  Freqs     <- seq(0.05, 0.95, 0.1)
  RepFragm  <- c(25, 25, 9, 9, 9, 13, 3, 3, 3, 0)
  RepFull   <- c(8, 3, 21, 16, 16, 5, 8, 5, 5, 13)
  FreqFragm <- unlist(lapply(1:length(Freqs), function(i) rep(Freqs[i], RepFragm[i])))
  FreqFull  <- unlist(lapply(1:length(Freqs), function(i) rep(Freqs[i], RepFull[i])))
}

cat("done!\n")


###################################################
#                                                 #
#  Explore fit for simulations                    #
#                                                 #
###################################################


cat("******  Exploring  grid    ***********\n\n")
# Results for distribution of full-length L1
ResultListFull_1000G_hist <- ExploreSelectionParameterGrid(FreqFull, 
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
                              SummaryType = SummaryType)

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
ResultListFragm_1000G_hist <- ExploreSelectionParameterGrid(FreqFragm, 
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
                               SummaryType = SummaryType)

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

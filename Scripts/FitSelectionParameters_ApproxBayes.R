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

# Path to ouput file
OutPath <- "D:/L1polymORF/Data/SelectionParameterFit_1000G_maxF0.6_Acc_Quant2.RData"

# Set population size
n = 10^4

# Set the number of replicates 
NrRep <- 10000

# Set the number of generations
NrGen <- 1000

# alpha and selection values for fine grid
aValsBasic_fineGrid = seq(51, 501, 10)
propS_fineGrid      = seq(0.70, 1, 0.02)

# Parameters of the gamma distribution of selection coefficients
Gshape = 90
GSscale = 1/100

# Maximum frequency (because high frequency insertions are in reference but not 
# 1000 genomes)
MaxF <- 0.6

# Probability vector for quantiles
PV <- seq(0, 0.40, 0.1)


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

# Load L1 catalog GenomicRanges
load("D:/L1polymORF/Data/L1CatalogGRanges.RData")
load('D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData')
cor.test(L1Catalogue$Allele_frequency_Num, L1Catalogue$ActivityNum)
sum(!is.na(L1Catalogue$Allele_frequency_Num))

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

# Count L1 per frequency bin
BinBreaks <- seq(0, 1, 0.2)
L1HistFragm <- hist(FreqFragm_1000G, breaks = BinBreaks, plot  = F)
L1HistFull <- hist(FreqFull_1000G,   breaks = BinBreaks, plot  = F)
par(mfrow = c(1, 1))
plot(L1HistFragm$mids, L1HistFull$density / L1HistFragm$density)
lines(c(0, 1), c(1, 1), lty = 2)


# Get the distance between catalog and 1000 Genome L1
DistCat2_1000G <- Dist2Closest(L1CatalogGR, L1_1000G_GRList_hg38$LiftedRanges)
DistCat2_1000G_hg19 <- Dist2Closest(L1CatalogGR_hg19, L1_1000G_GR_hg19)

# Get indices of 1000 Genome and catalog elements that match
idx1000G <- nearest(L1CatalogGR, L1_1000G_GRList_hg38$LiftedRanges)
L1CatalogMatch1000G <- L1CatalogL1Mapped[DistCat2_1000G < 100, ]
idx1000GMatchCat    <- idx1000G[DistCat2_1000G < 100]

# Result for distribution of catalog elements
blnCatFreq      <- !is.na(L1Catalogue$Allele_frequency_Num)
blnCatAct       <- L1Catalogue$ActivityNum > 0
blnCatFreq1000G <- is.na(L1CatalogMatch1000G$Allele_frequency)
blnCatAct1000G  <- L1CatalogMatch1000G$ActivityNum > 0
L1CatFreq <- c(L1Catalogue$Allele_frequency_Num[which(blnCatFreq & blnCatAct)],
               L1_1000G_reduced$Frequency[idx1000GMatchCat[blnCatFreq1000G & blnCatAct1000G]])
ks.test(L1Catalogue$Allele_frequency_Num[blnCatFreq],
        L1_1000G_reduced$Frequency[idx1000GMatchCat[blnCatFreq1000G]])
mean(L1Catalogue$Allele_frequency_Num[blnCatFreq])
mean(L1_1000G_reduced$Frequency[idx1000GMatchCat[blnCatFreq1000G]])
hist(L1_1000G_reduced$Frequency[idx1000GMatchCat[blnCatFreq1000G]])
mean(L1CatFreq)

cat("done!\n")

###################################################
#                                                 #
#  Functions to simulate allele frequencies       #
#                                                 #
###################################################

# Function to calculate Kolmogorov-Smirnov Statistic for the difference between
# simulated and observed frequencies.
DiffAlleleFreqKS <- function(ObservedFreq, Gshape, GSscale, n = 10^4, NrGen = 10^3,
                             NrRep = 10^4, NrSamples = 10^3, MaxFreq = 0.5){ 
  
  # Simulate allele frequencies
  AlleleFreq <- GenerateAlleleFreq(Gshape, GSscale, n = n, NrGen = NrGen,
                                   NrRep = NrRep, NrSamples = NrSamples)
  AlleleFreq[AlleleFreq>MaxFreq] <- 1
  
  # Calculate Kolmogorov-Smirnov statistic for the difference
  cat("Calculating Kolmogorov-Smirnov test\n\n")
  ks.test(AlleleFreq, ObservedFreq)$statistic
}

# Function to calculate differences in quantiles of simulated allele 
# frequencies.
QuantSimAlleleFreq <- function(Gshape, GSscale, n = 10^4, NrGen = 10^3,
                               NrRep = 10^4, NrSamples = 10^3, 
                               ProbV = seq(0.05, 0.95, 0.1), MaxFreq){ 
  
  # Simulate allele frequencies
  AlleleFreq <- GenerateAlleleFreq(Gshape, GSscale, n = n, NrGen = NrGen,
                                   NrRep = NrRep, NrSamples = NrSamples)
  AlleleFreq[AlleleFreq > MaxFreq] <- 1
  
  # Calculate qunatiles for simulated frequencies
  cat("Calculating quantiles\n\n")
  quantile(AlleleFreq, ProbV)
}

# # Function to explore a grid of alpha values and fitness values
# # Grid of parameter values
# ExploreGrid <- function(ObservedFreq, 
#                         aValsBasic = seq(1, 101, 10),
#                         proPs      = seq(0.5, 1.5, 0.1),
#                         PopSize  = 10^4,
#                         MaxFreq = 0.5,
#                         blnPlot = F){
#   
#   # Repeat proportion (fitness) values so that each alpha gets combined with
#   # the full range of fitness values
#   proPsRep   <- rep(proPs, length(aValsBasic))
#   aVals      <- rep(aValsBasic, each = length(proPs))
#   bVals      <- 1/aVals * proPsRep
#   
#   ObservedFreq[ObservedFreq > MaxFreq] <- 1
#   
#   # Evaluate difference according to Kolmogorov-Smirnov statistic 
#   DiffKS <- sapply(1:length(aVals), function(i){
#     DiffAlleleFreqKS(ObservedFreq, aVals[i], bVals[i], n = PopSize,
#                      MaxFreq = MaxFreq)
#   })
#   
#   # get indices of minimum difference per alpha value
#   idxMinDiffKS <- sapply(aValsBasic, function (x) {
#     idxA <- which(aVals == x)
#     idxMin <- which.min(DiffKS[idxA])
#     idxA[idxMin]})
#   
#   
# 
#   # Return values in a list
#   list(proPsRep  = proPsRep, aVals = aVals, bVals = bVals, DiffKS = DiffKS,
#        idxMinDiffKS = idxMinDiffKS)
# }

ExploreGrid <- function(ObservedFreq, 
                        aValsBasic = seq(1, 101, 10),
                        proPs      = seq(0.5, 1.5, 0.1),
                        PopSize  = 10^4,
                        MaxFreq = 0.5,
                        blnPlot = F,
                        SummaryType = c("none", "ks", "quant", "hist")){
  
  # Repeat proportion (fitness) values so that each alpha gets combined with
  # the full range of fitness values
  proPsRep   <- rep(proPs, length(aValsBasic))
  aVals      <- rep(aValsBasic, each = length(proPs))
  bVals      <- 1/aVals * proPsRep
  
  # Replace frequencies above maxium value by 1
  ObservedFreq[ObservedFreq > MaxFreq] <- 1
  
  # Determine the type of summary
  SummaryType <- SummaryType[1]
  
  # Loop over grid values and calculate summaries
  GridSummary <- switch(SummaryType,
    'none' = sapply(1:length(aVals), function(i){
                GenerateAlleleFreq(aVals[i], bVals[i],  n = n, NrGen = NrGen,
                         NrRep = NrRep, NrSamples = NrSamples)
             }),
    'ks'  = sapply(1:length(aVals), function(i){
              # Simulate allele frequencies
              AlleleFreq <- GenerateAlleleFreq(aVals[i], bVals[i],  n = n, NrGen = NrGen,
                                               NrRep = NrRep, NrSamples = NrSamples)
              AlleleFreq[AlleleFreq>MaxFreq] <- 1
      
              # Calculate Kolmogorov-Smirnov statistic for the difference
              cat("Calculating Kolmogorov-Smirnov test\n\n")
              ks.test(AlleleFreq, ObservedFreq)$statistic
            }),
    'quant' = sapply(1:length(aVals), function(i){
                 # Simulate allele frequencies
                 AlleleFreq <- GenerateAlleleFreq(Gshape, GSscale, n = n, NrGen = NrGen,
                                       NrRep = NrRep, NrSamples = NrSamples)
                 AlleleFreq[AlleleFreq > MaxFreq] <- 1
      
                # Calculate qunatiles for simulated frequencies
                cat("Calculating quantiles\n\n")
                quantile(AlleleFreq, ProbV)
            }),
    'hist' = sapply(1:length(aVals), function(i){
                # Simulate allele frequencies
                AlleleFreq <- GenerateAlleleFreq(Gshape, GSscale, n = n, NrGen = NrGen,
                                       NrRep = NrRep, NrSamples = NrSamples)
                AlleleFreq[AlleleFreq > MaxFreq] <- 1
      
                # Calculate qunatiles for simulated frequencies
                cat("Calculating histogram\n\n")
                hist(AlleleFreq, BreakV)$density
         })
    )
  

  # Return values in a list
  list(proPsRep  = proPsRep, aVals = aVals, bVals = bVals, DiffKS = DiffKS,
       idxMinDiffKS = idxMinDiffKS)
}

# Function to explore a grid of alpha values and fitness values
# and calculate qunatile differences
ExploreGrid_Quant <- function(ObservedFreq, 
                              aValsBasic = seq(1, 101, 10),
                              proPs      = seq(0.5, 1.5, 0.1),
                              PopSize  = 10^4,
                              blnPlot = F,
                              MaxFreq = 0.5,
                              ProbV = seq(0.05, 0.95, 0.1),
                              Epsilon = 0.1){
  
  # Repeat proportion (fitness) values so that each alpha gets combined with
  # the full range of fitness values
  proPsRep   <- rep(proPs, length(aValsBasic))
  aVals      <- rep(aValsBasic, each = length(proPs))
  bVals      <- 1/aVals * proPsRep
  
  # Evaluate difference according to Kolmogorov-Smirnov statistic 
  SimQuant <- sapply(1:length(aVals), function(i){
    QuantSimAlleleFreq(Gshape = aVals[i], GSscale = bVals[i], n = PopSize,
                       ProbV = ProbV, MaxFreq = MaxFreq)
  })
  ObsQuant      <- quantile(ObservedFreq, ProbV)
  DiffQuant     <- SimQuant - ObsQuant
  AbsDiffQuant  <- abs(DiffQuant)
  DiffQuantMean <- colMeans(AbsDiffQuant)
  DiffQuantMax  <- apply(AbsDiffQuant, 2, max)
  
  # Estimate intercept via regression
  blnE    <- DiffQuantMax < Epsilon
  XMat    <- t(DiffQuant)[blnE, ]
  if (sum(blnE) > 5){
    LMFit_a <- lm(aVals[blnE]    ~ XMat)
    LMFit_p <- lm(proPsRep[blnE] ~ XMat)
    aSample <- XMat %*% LMFit_a$coefficients[-1]
    pSample <- XMat %*% LMFit_p$coefficients[-1]
  } else {
    warning("Not enough values crossed threshold. No regression!\n")
    LMFit_a <- NA
    LMFit_p <- NA
    aSample <- NA
    pSample <- NA
    
  }
  
  # Return values in a list
  list(proPsRep  = proPsRep, aVals = aVals, bVals = bVals, 
       SimQuant = SimQuant,
       DiffQuant = DiffQuant, DiffQuantMean = DiffQuantMean,
       LMFit_a = LMFit_a, LMFit_p = LMFit_p, aSample = aSample,
       pSample = pSample)
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
# ResultList1Full_1000G <- ExploreGrid(FreqFull_1000G,
#                            aValsBasic = seq(1, 101, 10),
#                            proPs = seq(0.2, 1.5, 0.1))
# 
# # Result for distribution of catalog elements
# # blnCatFeq <- !is.na(L1Catalogue$Allele_frequency_Num)
# # ResultList1Cat <- ExploreGrid(L1Catalogue$Allele_frequency_Num[blnCatFeq],
# #                               aValsBasic = seq(1, 101, 10),
# #                               proPs = seq(0.2, 1.5, 0.1))
# 
# # Results for distribution of fragment L1
# ResultList1Fragm_1000G <- ExploreGrid(FreqFragm_1000G,
#                                aValsBasic = seq(1, 101, 10),
#                                proPs = seq(0.2, 1.5, 0.1))
# 
# 
# #######
# # Analyze finer grid
# #######

cat("******  Exploring fine grid    ***********\n\n")
# Results for distribution of full-length L1
ResultList2Full_1000G <- ExploreGrid_Quant(FreqFull_1000G,
                                           aValsBasic = aValsBasic_fineGrid,
                                           proPs = propS_fineGrid,
                                           MaxFreq = MaxF,
                                           ProbV = PV)

# Results for distribution of fragment L1
ResultList2Fragm_1000G <- ExploreGrid_Quant(FreqFragm_1000G,
                                            aValsBasic = aValsBasic_fineGrid,
                                            proPs = propS_fineGrid,
                                            MaxFreq = MaxF,
                                            ProbV = PV)

#########
# Save results
#########

save.image(OutPath)

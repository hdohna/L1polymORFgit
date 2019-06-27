# The script below estimates the effect of L1 width on selection coefficients
# from 1000 genome data for different populations 
# 

# The script below estimates selection coefficients of L1 from the 
# 1000 genome data 

##########################################
#                                        #
#     Load packages                      #
#                                        #
##########################################

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(GenomicRanges)

##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

# Specify file paths
DataPath        <- 'D:/L1polymORF/Data/'
G1000SamplePath <- 'D:/L1polymORF/Data/1000GenomeSampleInfo.txt'
L1GRPath        <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
L1RefRangePath  <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
ChrLPath        <- 'D:/L1polymORF/Data/ChromLengthsHg19.Rdata'
L1RefPath       <- 'D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv'
SelectResultOutPath <- "D:/L1polymORF/Data/L1SelectionResults_byPop.RData"

# False discovery rate for selected L1
FDR <- 0.1

# Specify range width for DNAse analysis
RangeWidth <- 10^6

# Human effective population size
PopSize <- 10^5

# Minimum length for a full L1
MinLengthFullL1 <- 6000

# Sample size for ME insertion calls
MEInsSamplesize <- 2504

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("\n\nLoading and processing data ...")

# Load previously generated objects
load(L1GRPath)
load(ChrLPath)
load(L1RefRangePath)

# Read table with sample info and get IDs belonging to particular populations
SampleInfo <- read.table(G1000SamplePath, header = T, as.is = T)

# Calculate population frequency
PopFreq <- sapply(unique(SampleInfo$super_pop), function(x){
  blnPop <- SampleInfo$super_pop == x
  rowSums(L1_1000G[,SampleInfo$sample[blnPop]])
})

# Create dataframe with 
PredictMat <- data.frame(blnFull = L1_1000G$InsLength >= 6000,
                        L1width = L1_1000G$InsLength)
PredictMat <- cbind(PredictMat, PopFreq)

cat("done!\n")

###################################################
#                                                 #
#   Fit effect of insertion length on selection   #
#                                                 #
###################################################

# Create a matrix of predictor variables (L1 start and boolean variable for)
blnNA      <- sapply(1:nrow(PredictMat), function(x) any(is.na(PredictMat[x,])))
PredictMat <- PredictMat[!blnNA, ]

cat("\n********   Estimating effect of insertion length: Europeans    **********\n")
blnEur      <- PredictMat$EUR > 0
ModelFitEUR <- FitSelectionModels(PredictMat[blnEur,1:3],  
                                Freqs = round(PredictMat$EUR[blnEur], 0), 
                                Counts = rep(1, sum(blnEur)), 
                                PopSize = PopSize, 
                                SampleSize = sum(SampleInfo$super_pop == "EUR"),
                                blnIns = rep(T, sum(blnEur)), 
                                LogRegCoeff = c(-4.706691, 9.736755),
                                DetectProb = rep(0.8, sum(blnEur)),
                                aBorder = 0.002, 
                                bBorder = 10^(-2), 
                                cBorder = 10^(-5))

cat("\n********   Estimating effect of insertion length: East Asians   **********\n")
blnEAS      <- PredictMat$EAS > 0
ModelFitEAS <- FitSelectionModels(PredictMat[blnEAS,1:3],  
                                  Freqs = round(PredictMat$EAS[blnEAS], 0), 
                                  Counts = rep(1, sum(blnEAS)), 
                                  PopSize = PopSize, 
                                  SampleSize = sum(SampleInfo$super_pop == "EAS"),
                                  blnIns = rep(T, sum(blnEAS)), 
                                  LogRegCoeff = c(-4.706691, 9.736755),
                                  DetectProb = rep(0.8, sum(blnEAS)),
                                  aBorder = 0.003, 
                                  bBorder = 10^(-2), 
                                  cBorder = 10^(-5))

cat("\n********   Estimating effect of insertion length: Americans   **********\n")
blnAMR      <- PredictMat$AMR > 0
ModelFitAMR <- FitSelectionModels(PredictMat[blnAMR,1:3],  
                                  Freqs = round(PredictMat$AMR[blnAMR], 0), 
                                  Counts = rep(1, sum(blnAMR)), 
                                  PopSize = PopSize, 
                                  SampleSize = sum(SampleInfo$super_pop == "AMR"),
                                  blnIns = rep(T, sum(blnAMR)), 
                                  LogRegCoeff = c(-4.706691, 9.736755),
                                  DetectProb = rep(0.8, sum(blnAMR)),
                                  aBorder = 0.003, 
                                  bBorder = 10^(-2), 
                                  cBorder = 10^(-5))

cat("\n********   Estimating effect of insertion length: South Asians   **********\n")
blnSAS      <- PredictMat$SAS > 0
ModelFitSAS <- FitSelectionModels(PredictMat[blnSAS,1:3],  
                                  Freqs = round(PredictMat$SAS[blnSAS], 0), 
                                  Counts = rep(1, sum(blnSAS)), 
                                  PopSize = PopSize, 
                                  SampleSize = sum(SampleInfo$super_pop == "SAS"),
                                  blnIns = rep(T, sum(blnSAS)), 
                                  LogRegCoeff = c(-4.706691, 9.736755),
                                  DetectProb = rep(0.8, sum(blnSAS)),
                                  aBorder = 0.003, 
                                  bBorder = 10^(-2), 
                                  cBorder = 10^(-5))


cat("\n********   Estimating effect of insertion length: Africans   **********\n")
blnAFR      <- PredictMat$AFR > 0
ModelFitAFR <- FitSelectionModels(PredictMat[blnAFR,1:3],  
                                  Freqs = round(PredictMat$AFR[blnAFR], 0), 
                                  Counts = rep(1, sum(blnAFR)), 
                                  PopSize = PopSize, 
                                  SampleSize = sum(SampleInfo$super_pop == "AFR"),
                                  blnIns = rep(T, sum(blnAFR)), 
                                  LogRegCoeff = c(-4.706691, 9.736755),
                                  DetectProb = rep(0.8, sum(blnAFR)),
                                  aBorder = 0.003, 
                                  bBorder = 10^(-2), 
                                  cBorder = 10^(-5))


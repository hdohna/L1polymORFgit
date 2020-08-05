# The script below estimates selection coefficients of L1 from the 
# 1000 genome data 

##########################################
#                                        #
#     Load packages                      #
#                                        #
##########################################

# Load packages
library(GenomicRanges)
library(pracma)

# Source start script
source('D:/OneDrive - American University of Beirut/L1polymORFgit/Scripts/_Start_L1polymORF.R')


##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

# Load previously generated objects
#load('D:/OneDrive - American University of Beirut/L1polymORF/Data/SingletonAnalysis_unphased.RData')

# Load simulation results
load("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1Simulated_AdditionalInfo_MELT.RData")

# Specify file paths
DataPath            <- 'D:/OneDrive - American University of Beirut/L1polymORF/Data/'
L1GRPath            <-  "D:/OneDrive - American University of Beirut/L1polymORF/Data/GRanges_L1_1000Genomes.RData"
MeltInsPath         <- "D:/OneDrive - American University of Beirut/L1polymORF/Data/nstd144.GRCh37.variant_call.vcf"
MeltDelPath         <- "D:/OneDrive - American University of Beirut/L1polymORF/Data/DEL.final_comp.vcf"
ChrLPath            <- 'D:/OneDrive - American University of Beirut/L1polymORF/Data/ChromLengthsHg19.Rdata'
L1RefPath           <- 'D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HS_repeat_table_Hg19.csv'
L1RefRangePath      <- 'D:/OneDrive - American University of Beirut/L1polymORF/Data/L1RefRanges_hg19.Rdata'
RegrOutputPath      <- "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1RegressionResults.RData"
SelectResultOutPath <- "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1SelectionResults_MELT_GroupwithSim.RData"

# False discovery rate for selected L1
FDR <- 0.1

# Specify range width for DNAse analysis
RangeWidth <- 10^6

# Human effective population size
PopSize <- 10^5

# Minimum length for a full L1
MinLengthFullL1 <- 6000

# Sample size for ME insertion calls
MEInsSamplesize <- 2453

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("\n\nLoading and processing data ...")

# Source start script again
source('D:/OneDrive - American University of Beirut/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Read in vcf file with MELT insertion calls
MEInsCall <- read.table(MeltInsPath, 
                        as.is = T,
                        col.names = c("Chrom", "Pos", "ID", "Alt", "Type", "V6", 
                                      "V7", "Info"))
MEInsCall <- MEInsCall[MEInsCall$Type == "<INS:ME:LINE1>",]
grep("L1Ta1", MEInsCall$Info)

# Extract allele frequency from info column
GetAF <- function(x){
  xSplit <- strsplit(x, ";")[[1]]
  AFch   <- strsplit(xSplit[length(xSplit)], "=")[[1]][2]
  as.numeric(AFch)
}
GetLength <- function(x){
  xSplit <- strsplit(x, ";")[[1]]
  LengthCh   <- strsplit(xSplit[grep("SVLEN=", xSplit)], "=")[[1]][2]
  as.numeric(LengthCh)
}

# Add columns necessary for analysis 
MEInsCall$AF <- sapply(MEInsCall$Info, GetAF)
MEInsCall <- MEInsCall[!is.na(MEInsCall$AF), ]
MEInsCall$L1width <- sapply(MEInsCall$Info, GetLength)
MEInsCall$SampleSize <- 1/min(MEInsCall$AF) 
# MEInsCall$SampleSize <- 2 * MEInsSamplesize
MEInsCall$Freq <- ceiling(MEInsCall$SampleSize * MEInsCall$AF) # TODO: Figure out why not integers!
MEInsCall$blnFull <- MEInsCall$L1width >= MinLengthFullL1
length(unique(MEInsCall$Freq))
# Create GRanges object for MEInsCall
MEInsCall$ChromName <- paste("chr", MEInsCall$Chrom, sep = "")
MEIns_GR <- makeGRangesFromDataFrame(df = MEInsCall,
                                     seqnames.field = "ChromName",
                                     start.field = "Pos",
                                     end.field = "Pos")

# Read in vcf file with MELT deletion calls
MEDelCall <- ReadVCF(MeltDelPath)
MEDelCall$chromosome <- paste("chr", MEDelCall$X.CHROM, sep = "")
MEDel_GR  <- makeGRangesFromDataFrame(df = MEDelCall,
                                     start.field = "POS",
                                     end.field = "POS")
colnames(MEDelCall)
MEDelCall$INFO[1:5]

# function to get numeric genotype
GetNumericGenotype <- function(x){
  Split1 <- strsplit(x, ":")[[1]][1]
  Split2 <- strsplit(Split1, "/")[[1]]
  sum(as.numeric(Split2))
}

# Get numeric genotype of all reference L1 deletions
GTCols <- grep("L1Filtered", colnames(MEDelCall))
L1RefNumGen <- 2 - sapply(GTCols, function(x){
  sapply(1:nrow(MEDelCall), function(y) GetNumericGenotype(MEDelCall[y,x]))
})

# Add columns for frequency and sample size
MEDelCall$Freq       <- rowSums(L1RefNumGen, na.rm = T)
MEDelCall$SampleSize <- apply(L1RefNumGen, 1, function(x) 2*sum(!is.na(x)))

load(L1GRPath)
load(ChrLPath)
load(L1RefRangePath)
load(RegrOutputPath)
load("D:/OneDrive - American University of Beirut/L1polymORF/Data/DelVsL1Length.RData")

# Load table with L1 reference data
L1RefTable <- read.csv(L1RefPath, as.is = T)

# Create genomic ranges of reference L1 with 100 bp added on each side
L1NeighborRanges <- GRanges(seqnames = seqnames(L1GRanges), 
                            IRanges(start = start(L1GRanges) - 100,
                                    end   = end(L1GRanges) + 100))

# Create a data frame of reference L1
RefL1Data <- data.frame(L1width = width(L1GRanges),
                        Freq = 30, SampleSize = 30)
OL_MEDelRefL1 <- findOverlaps(L1NeighborRanges, MEDel_GR)
RefL1Data$Freq[OL_MEDelRefL1@from] <- MEDelCall$Freq[OL_MEDelRefL1@to]
RefL1Data$SampleSize[OL_MEDelRefL1@from] <- MEDelCall$SampleSize[OL_MEDelRefL1@to]
RefL1Data$blnFull <- RefL1Data$L1width >= MinLengthFullL1
RefL1Data$Info    <- NA

# Number of L1 that are fixed at proportion 1
N1 <- length(L1GRanges) - length(OL_MEDelRefL1@from)

RefL1Data <- RefL1Data[OL_MEDelRefL1@from, ]
L1GRanges <- L1GRanges[OL_MEDelRefL1@from]

# Put data of non-reference L1 (insertions) and reference L1 (deletions) 
# together
L1TotData <- rbind(MEInsCall[ ,c("L1width", "Freq", "SampleSize", "blnFull",
                                 "Info")],
                   RefL1Data)
L1TotData$blnIns <- c(rep(T, nrow(MEInsCall)), rep(F, nrow(RefL1Data)))
L1TotData$L1Freq <- NA
L1TotData$L1Freq[L1TotData$blnIns] <- L1TotData$Freq[L1TotData$blnIns] / 
  L1TotData$SampleSize[L1TotData$blnIns]
L1TotData$L1Freq[!L1TotData$blnIns] <- 1 - L1TotData$Freq[!L1TotData$blnIns] / 
  L1TotData$SampleSize[!L1TotData$blnIns]
L1TotData$DetectProb <- 0.85
L1TotData$DetectProb[L1TotData$blnIns] <- 0.9

# Perform logistic regression for the probability of reference L1 as function
# of L1 frequency
L1TotData$blnRef <- !L1TotData$blnIns
LogRegL1Ref <- glm(blnRef ~ L1Freq, family = binomial, data = L1TotData)
LogRegL1Ref$coefficients
L1TotData$L1Freq
# Combine genomic ranges
L1TotGR <- c(MEIns_GR, L1GRanges)

# Number of L1 that are not fixed
Nnf <- nrow(L1TotData)

# Add sampled L1 width
L1DetectAgg_withL1 <- L1DetectAgg_Group[
  which(L1DetectAgg_Group$L1widthTrue > 0 & L1DetectAgg_Group$L1widthEst > 0 ),]
blnL1widthNA <- is.na(L1TotData$L1width)
SampleTrueL1Width(
  SimL1widthTrue = L1DetectAgg_withL1$L1widthTrue, 
  SimL1widthEst  = L1DetectAgg_withL1$L1widthEst,
  SimL1Freq = L1DetectAgg_withL1$EstFreq,
  EstL1width = L1TotData$L1width[!blnL1widthNA],
  ObsL1Freq = L1TotData$Freq[!blnL1widthNA],
  L1widthBreaks = seq(0, 6500, 500), PlotPath = "D:/OneDrive - American University of Beirut/L1polymORF/Figures/PL1TrueWidthVsFreq.pdf")
cat("done!\n")

# Test for frequency difference between L1 labelled Ta and not Ta
t.test(L1TotData$Freq[grep("L1Ta", L1TotData$Info)],
       L1TotData$Freq[-grep("L1Ta", L1TotData$Info)])
wilcox.test(L1TotData$Freq[grep("L1Ta", L1TotData$Info)],
       L1TotData$Freq[-grep("L1Ta", L1TotData$Info)])


###################################################
#                                                 #
#   Fit effect of insertion length on selection   #
#                                                 #
###################################################

# Create a matrix of predictor variables (L1 start and boolean variable for)
PredictMat <- L1TotData[, c("blnFull", "L1width", "Freq")]
blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))
sum(blnNA)

# AlleleFreqLogLik_4Par_pracma(
#   Freqs = round(L1TotData$Freq[!blnNA], 0), 
#   Counts = rep(1, sum(!blnNA)), 
#   Predict = PredictMat[!blnNA, 1:3],
#   a =  ModelFit_pracma$ML_abc$`par`[1], b = ModelFit_pracma$ML_abc$`par`[2], 
#   c = ModelFit_pracma$ML_abc$`par`[3], d = 0, N = PopSize, 
#   SampleSize = L1TotData$SampleSize[!blnNA],
#   blnIns = L1TotData$blnIns[!blnNA], 
#   LogRegCoeff = LogRegL1Ref$coefficients,
#   DetectProb = L1TotData$DetectProb[!blnNA])

cat("\n********   Estimating effect of insertion length: true length   **********\n")
# ModelFit_romberg <- FitSelectionModels_romberg(PredictMat[!blnNA, 1:3],  
#                                              Freqs = round(L1TotData$Freq[!blnNA], 0), 
#                                              Counts = rep(1, sum(!blnNA)), 
#                                              PopSize = PopSize, 
#                                              SampleSize = L1TotData$SampleSize[!blnNA],
#                                              blnIns = L1TotData$blnIns[!blnNA], 
#                                              LogRegCoeff = LogRegL1Ref$coefficients,
#                                              DetectProb = L1TotData$DetectProb[!blnNA],
#                                              aBorder = 0.003, 
#                                              bBorder = 10^(-2), 
#                                              cBorder = 10^(-5))
# ModelFit_quadgr <- FitSelectionModels_quadgr(PredictMat[!blnNA, 1:3],  
#                                                Freqs = round(L1TotData$Freq[!blnNA], 0), 
#                                                Counts = rep(1, sum(!blnNA)), 
#                                                PopSize = PopSize, 
#                                                SampleSize = L1TotData$SampleSize[!blnNA],
#                                                blnIns = L1TotData$blnIns[!blnNA], 
#                                                LogRegCoeff = LogRegL1Ref$coefficients,
#                                                DetectProb = L1TotData$DetectProb[!blnNA],
#                                                aBorder = 0.003, 
#                                                bBorder = 10^(-2), 
#                                                cBorder = 10^(-5))

cat("\n********   Estimating effect of insertion length: true length   **********\n")
ModelFit_pracma <- FitSelectionModels_pracma(PredictMat[!blnNA, 1:3],  
                                Freqs = round(L1TotData$Freq[!blnNA], 0), 
                                Counts = rep(1, sum(!blnNA)), 
                                PopSize = PopSize, 
                                SampleSize = L1TotData$SampleSize[!blnNA],
                                blnIns = L1TotData$blnIns[!blnNA], 
                                LogRegCoeff = LogRegL1Ref$coefficients,
                                DetectProb = L1TotData$DetectProb[!blnNA],
                                aBorder = 0.003, 
                                bBorder = 10^(-2), 
                                cBorder = 10^(-5))

cat("\n********   Estimating effect of insertion length: true length   **********\n")
# ModelFit1 <- FitSelectionModels(PredictMat[!blnNA, 1:3],  
#                                 Freqs = round(L1TotData$Freq[!blnNA], 0), 
#                                 Counts = rep(1, sum(!blnNA)), 
#                                 PopSize = PopSize, 
#                                 SampleSize = L1TotData$SampleSize[!blnNA],
#                                 blnIns = L1TotData$blnIns[!blnNA], 
#                                 LogRegCoeff = LogRegL1Ref$coefficients,
#                                 DetectProb = L1TotData$DetectProb[!blnNA],
#                                 aBorder = 0.003, 
#                                 bBorder = 10^(-2), 
#                                 cBorder = 10^(-5))

# # Fit model using only L1 Ta-1 
# PredictMat <- L1TotData[grep("L1Ta1", L1TotData$Info), c("blnFull", "L1width", "Freq", "Info")]
# blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))
# 
# cat("\n********   Estimating effect of insertion length: true length for TA1   **********\n")
# ModelFit_Ta1 <- FitSelectionModels(PredictMat[!blnNA, 1:3],  
#                                 Freqs = round(L1TotData$Freq[!blnNA], 0), 
#                                 Counts = rep(1, sum(!blnNA)), 
#                                 PopSize = PopSize, 
#                                 SampleSize = L1TotData$SampleSize[!blnNA],
#                                 blnIns = L1TotData$blnIns[!blnNA], 
#                                 LogRegCoeff = LogRegL1Ref$coefficients,
#                                 DetectProb = L1TotData$DetectProb[!blnNA],
#                                 aBorder = 0.003, 
#                                 bBorder = 10^(-2), 
#                                 cBorder = 10^(-5))
# 
# # Fit model using only L1 Ta-1 
# idxNoTa1 <- setdiff(grep("L1Ta", L1TotData$Info), grep("L1Ta1", L1TotData$Info))
# PredictMat <- L1TotData[idxNoTa1, c("blnFull", "L1width", "Freq", "Info")]
# blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))
# 
# cat("\n********   Estimating effect of insertion length: true length for TA1   **********\n")
# ModelFit_noTa1 <- FitSelectionModels(PredictMat[!blnNA, 1:3],  
#                                    Freqs = round(L1TotData$Freq[!blnNA], 0), 
#                                    Counts = rep(1, sum(!blnNA)), 
#                                    PopSize = PopSize, 
#                                    SampleSize = L1TotData$SampleSize[!blnNA],
#                                    blnIns = L1TotData$blnIns[!blnNA], 
#                                    LogRegCoeff = LogRegL1Ref$coefficients,
#                                    DetectProb = L1TotData$DetectProb[!blnNA],
#                                    aBorder = 0.003, 
#                                    bBorder = 10^(-2), 
#                                    cBorder = 10^(-5))
# # Fit model using only L1 Ta 
# PredictMat <- L1TotData[grep("L1Ta", L1TotData$Info), c("blnFull", "L1width", "Freq", "Info")]
# blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))
# 
# cat("\n********   Estimating effect of insertion length: true length for Ta   **********\n")
# ModelFit_Ta <- FitSelectionModels(PredictMat[!blnNA, 1:3],  
#                                    Freqs = round(L1TotData$Freq[!blnNA], 0), 
#                                    Counts = rep(1, sum(!blnNA)), 
#                                    PopSize = PopSize, 
#                                    SampleSize = L1TotData$SampleSize[!blnNA],
#                                    blnIns = L1TotData$blnIns[!blnNA], 
#                                    LogRegCoeff = LogRegL1Ref$coefficients,
#                                    DetectProb = L1TotData$DetectProb[!blnNA],
#                                    aBorder = 0.003, 
#                                    bBorder = 10^(-2), 
#                                    cBorder = 10^(-5))
# 
# # Fit model using only L1 non-Ta 
# PredictMat <- L1TotData[-grep("L1Ta", L1TotData$Info), c("blnFull", "L1width", "Freq", "Info")]
# blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))
# 
# cat("\n********   Estimating effect of insertion length: true length for Ta   **********\n")
# ModelFit_nonTa <- FitSelectionModels(PredictMat[!blnNA, 1:3],  
#                                   Freqs = round(L1TotData$Freq[!blnNA], 0), 
#                                   Counts = rep(1, sum(!blnNA)), 
#                                   PopSize = PopSize, 
#                                   SampleSize = L1TotData$SampleSize[!blnNA],
#                                   blnIns = L1TotData$blnIns[!blnNA], 
#                                   LogRegCoeff = LogRegL1Ref$coefficients,
#                                   DetectProb = L1TotData$DetectProb[!blnNA],
#                                   aBorder = 0.003, 
#                                   bBorder = 10^(-2), 
#                                   cBorder = 10^(-5))

# # Fit model including  an effect of full-length L1 that differs between Ta and non-Ta
# blnTa <- 1:nrow(L1TotData) %in% grep("L1Ta", L1TotData$Info)
# L1TotData$blnFull_Ta    <- L1TotData$blnFull & blnTa
# L1TotData$blnFull_nonTa <- L1TotData$blnFull & (!blnTa)
# L1TotData$blnFull_Ta    <- L1TotData$blnFull 
# L1TotData$blnFull_nonTa <- 0
# PredictMat <- L1TotData[, c("blnFull_Ta", "blnFull_nonTa", "L1width", "Freq")]
# blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))
# sum(blnNA)
# 
# cat("Estimate effect of L1 width and blnFull, spearated by Ta and non-Ta on selection ...\n")
# ML_L1WidthFullTa_nonTa <- constrOptim(theta = c(a = ModelFit_pracma$ML_abc$`par`[1], 
#                                 b = ModelFit_pracma$ML_abc$`par`[2], 
#                                 c = ModelFit_pracma$ML_abc$`par`[2],
#                                 d = ModelFit_pracma$ML_abc$`par`[3]),
#                                 
#                                 f = function(x) -AlleleFreqLogLik_4Par_pracma(
#                                   Freqs = round(L1TotData$Freq[!blnNA], 0), 
#                                   Counts = rep(1, sum(!blnNA)), 
#                                   Predict = PredictMat[!blnNA, 1:3],
#                                   a = x[1], b = x[2], c = x[3], d = x[4], N = PopSize, 
#                                   SampleSize = L1TotData$SampleSize[!blnNA],
#                                   blnIns = L1TotData$blnIns[!blnNA], 
#                                   LogRegCoeff = LogRegL1Ref$coefficients,
#                                   DetectProb = L1TotData$DetectProb[!blnNA]),
#                                 
#                       grad = NULL,
#                       ui = rbind(c(1, 0, 0, 0),  c(0, 1, 0, 0),  c(0, 0, 1, 0), c(0, 0, 0, 1),
#                                  c(-1, 0, 0, 0),  c(0, -1, 0, 0),  c(0, 0,- 1, 0), c(0, 0, 0, -1)),
#                       ci = c(a = -0.003, b1 = -10^(-2), b2 = -10^(-2), c = -10^(-5), 
#                              a = -0.003, b1 = -10^(-2), b2 = -10^(-2), c = -10^(-5)),
#                       method = "Nelder-Mead")
# # Function to extract AIC from optim results
# cat("done!\n")
# # Function to extract AIC from optim results
# GetAIC <- function(OptimResults){
#   round(2 * (length(OptimResults$par) + OptimResults$value), 2)
# }
# GetAIC(ML_L1WidthFullTa_nonTa)
# 
# # Fit model including an intercept differs between Ta and non-Ta
# L1TotData$blnTa <- 1:nrow(L1TotData) %in% grep("L1Ta", L1TotData$Info)
# PredictMat <- L1TotData[, c("blnFull", "blnTa", "L1width", "Freq")]
# blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))
# sum(blnNA)
# 
# cat("Estimate effect of L1 width, blnFull, spearated by Ta and non-Ta on selection ...\n")
# ML_L1WidthFull_InterceptTa_nonTa <- constrOptim(theta = c(a = ModelFit_pracma$ML_abc$`par`[1], 
#                                                 b = ModelFit_pracma$ML_abc$`par`[2], 
#                                                 c = 0,
#                                                 d = ModelFit_pracma$ML_abc$`par`[3]),
#                                       
#                                       f = function(x) -AlleleFreqLogLik_4Par_pracma(
#                                         Freqs = round(L1TotData$Freq[!blnNA], 0), 
#                                         Counts = rep(1, sum(!blnNA)), 
#                                         Predict = PredictMat[!blnNA, 1:3],
#                                         a = x[1], b = x[2], c = x[3], d = x[4], N = PopSize, 
#                                         SampleSize = L1TotData$SampleSize[!blnNA],
#                                         blnIns = L1TotData$blnIns[!blnNA], 
#                                         LogRegCoeff = LogRegL1Ref$coefficients,
#                                         DetectProb = L1TotData$DetectProb[!blnNA]),
#                                       
#                                       grad = NULL,
#                                       ui = rbind(c(1, 0, 0, 0),  c(0, 1, 0, 0),  c(0, 0, 1, 0), c(0, 0, 0, 1),
#                                                  c(-1, 0, 0, 0),  c(0, -1, 0, 0),  c(0, 0,- 1, 0), c(0, 0, 0, -1)),
#                                       ci = c(a = -0.003, b1 = -10^(-2), b2 = -10^(-2), c = -10^(-5), 
#                                              a = -0.003, b1 = -10^(-2), b2 = -10^(-2), c = -10^(-5)),
#                                       method = "Nelder-Mead")
# # Function to extract AIC from optim results
# cat("done!\n")
# # Function to extract AIC from optim results
# GetAIC <- function(OptimResults){
#   round(2 * (length(OptimResults$par) + OptimResults$value), 2)
# }
# GetAIC(ML_L1WidthFull_InterceptTa_nonTa)
# 
# # Fit model where relationship between L1 length and selection differs between Ta and non-Ta
# L1TotData$blnTa    <- 1:nrow(L1TotData) %in% grep("L1Ta", L1TotData$Info)
# L1TotData$blnTaL1w <- L1TotData$blnTa * L1TotData$L1width
# PredictMat <- L1TotData[, c("blnFull", "blnTaL1w", "L1width", "Freq")]
# blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))
# sum(blnNA)
# 
# cat("Estimate effect of L1 width, blnFull, spearated by Ta and non-Ta on selection ...\n")
# ML_L1WidthFull_SlopeTa_nonTa <- constrOptim(theta = c(a = ModelFit_pracma$ML_abc$`par`[1], 
#                                                           b = ModelFit_pracma$ML_abc$`par`[2], 
#                                                           c = 0,
#                                                           d = ModelFit_pracma$ML_abc$`par`[3]),
#                                                 
#                                                 f = function(x) -AlleleFreqLogLik_4Par_pracma(
#                                                   Freqs = round(L1TotData$Freq[!blnNA], 0), 
#                                                   Counts = rep(1, sum(!blnNA)), 
#                                                   Predict = PredictMat[!blnNA, 1:3],
#                                                   a = x[1], b = x[2], c = x[3], d = x[4], N = PopSize, 
#                                                   SampleSize = L1TotData$SampleSize[!blnNA],
#                                                   blnIns = L1TotData$blnIns[!blnNA], 
#                                                   LogRegCoeff = LogRegL1Ref$coefficients,
#                                                   DetectProb = L1TotData$DetectProb[!blnNA]),
#                                                 
#                                                 grad = NULL,
#                                                 ui = rbind(c(1, 0, 0, 0),  c(0, 1, 0, 0),  c(0, 0, 1, 0), c(0, 0, 0, 1),
#                                                            c(-1, 0, 0, 0),  c(0, -1, 0, 0),  c(0, 0,- 1, 0), c(0, 0, 0, -1)),
#                                                 ci = c(a = -0.003, b1 = -10^(-2), b2 = -10^(-2), c = -10^(-5), 
#                                                        a = -0.003, b1 = -10^(-2), b2 = -10^(-2), c = -10^(-5)),
#                                                 method = "Nelder-Mead")
# # Function to extract AIC from optim results
# cat("done!\n")
# # Function to extract AIC from optim results
# GetAIC <- function(OptimResults){
#   round(2 * (length(OptimResults$par) + OptimResults$value), 2)
# }
# GetAIC(ML_L1WidthFull_SlopeTa_nonTa)
# 
# # Fit model including  an effect of full-length L1 that differs between Ta1 and non-Ta1
# blnTa1 <- 1:nrow(L1TotData) %in% grep("L1Ta1", L1TotData$Info)
# L1TotData$blnFull_Ta1    <- L1TotData$blnFull & blnTa1
# L1TotData$blnFull_nonTa1 <- L1TotData$blnFull & (!blnTa1)
# L1TotData$blnFull_Ta1    <- L1TotData$blnFull 
# L1TotData$blnFull_nonTa1 <- 0
# PredictMat <- L1TotData[, c("blnFull_Ta1", "blnFull_nonTa1", "L1width", "Freq")]
# blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))
# sum(blnNA)
# 
# cat("Estimate effect of L1 width and blnFull, spearated by Ta1 and non-Ta1 on selection ...\n")
# ML_L1WidthFullTa1_nonTa1 <- constrOptim(theta = c(a = ModelFit_pracma$ML_abc$`par`[1], 
#                                                 b = ModelFit_pracma$ML_abc$`par`[2], 
#                                                 c = ModelFit_pracma$ML_abc$`par`[2],
#                                                 d = ModelFit_pracma$ML_abc$`par`[3]),
#                                       
#                                       f = function(x) -AlleleFreqLogLik_4Par_pracma(
#                                         Freqs = round(L1TotData$Freq[!blnNA], 0), 
#                                         Counts = rep(1, sum(!blnNA)), 
#                                         Predict = PredictMat[!blnNA, 1:3],
#                                         a = x[1], b = x[2], c = x[3], d = x[4], N = PopSize, 
#                                         SampleSize = L1TotData$SampleSize[!blnNA],
#                                         blnIns = L1TotData$blnIns[!blnNA], 
#                                         LogRegCoeff = LogRegL1Ref$coefficients,
#                                         DetectProb = L1TotData$DetectProb[!blnNA]),
#                                       
#                                       grad = NULL,
#                                       ui = rbind(c(1, 0, 0, 0),  c(0, 1, 0, 0),  c(0, 0, 1, 0), c(0, 0, 0, 1),
#                                                  c(-1, 0, 0, 0),  c(0, -1, 0, 0),  c(0, 0,- 1, 0), c(0, 0, 0, -1)),
#                                       ci = c(a = -0.003, b1 = -10^(-2), b2 = -10^(-2), c = -10^(-5), 
#                                              a = -0.003, b1 = -10^(-2), b2 = -10^(-2), c = -10^(-5)),
#                                       method = "Nelder-Mead")
# # Function to extract AIC from optim results
# cat("done!\n")
# 
# # Get AIC
# GetAIC(ML_L1WidthFullTa1_nonTa1) - GetAIC(ML_L1WidthFull_SlopeTa_nonTa)
# GetAIC(ML_L1WidthFull_InterceptTa_nonTa) - GetAIC(ML_L1WidthFull_SlopeTa_nonTa)
# as.numeric(as.character(ModelFit_pracma$AICTab$AIC)) - GetAIC(ML_L1WidthFull_SlopeTa_nonTa)
# 
# # # Fit model including  an effect of full-length L1 that differs between Ta1 and non-Ta1
# # blnTa1 <- 1:nrow(L1TotData) %in% grep("L1Ta1", L1TotData$Info)
# # L1TotData$blnFull_Ta1 <- L1TotData$blnFull & blnTa1
# # L1TotData$blnFull_nonTa1 <- L1TotData$blnFull & (!blnTa1)
# # PredictMat <- L1TotData[, c("blnFull_Ta1", "blnFull_nonTa1", "L1width", "Freq")]
# # blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))
# # 
# # cat("Estimate effect of L1 width and blnFull, spearated by Ta1 and non-Ta1 on selection ...\n")
# # 
# # ML_L1WidthFullTa1_nonTa1 <- constrOptim(theta = c(a = ModelFit1$ML_abc$`par`[1], 
# #                                                 b = ModelFit1$ML_abc$`par`[2], 
# #                                                 c = ModelFit1$ML_abc$`par`[2],
# #                                                 d = ModelFit1$ML_abc$`par`[3]),
# #                                       
# #                                       f = function(x) -AlleleFreqLogLik_4Par(
# #                                         Freqs = round(L1TotData$Freq[!blnNA], 0), 
# #                                         Counts = rep(1, sum(!blnNA)), 
# #                                         Predict = PredictMat[!blnNA, 1:3],
# #                                         a = x[1], b = x[2], c = x[3], d = x[4], N = PopSize, 
# #                                         SampleSize = L1TotData$SampleSize[!blnNA],
# #                                         blnIns = L1TotData$blnIns[!blnNA], 
# #                                         LogRegCoeff = LogRegL1Ref$coefficients,
# #                                         DetectProb = L1TotData$DetectProb[!blnNA]),
# #                                       
# #                                       grad = NULL,
# #                                       ui = rbind(c(1, 0, 0, 0),  c(0, 1, 0, 0),  c(0, 0, 1, 0), c(0, 0, 0, 1),
# #                                                  c(-1, 0, 0, 0),  c(0, -1, 0, 0),  c(0, 0,- 1, 0), c(0, 0, 0, -1)),
# #                                       ci = c(a = -0.003, b1 = -10^(-2), b2 = -10^(-2), c = -10^(-5), 
# #                                              a = -0.003, b1 = -10^(-2), b2 = -10^(-2), c = -10^(-5)),
# #                                       method = "Nelder-Mead")
# # GetAIC(ML_L1WidthFullTa1_nonTa1)
# # ML_L1WidthFullTa1_nonTa1$par
# # ML_L1WidthFullTa_nonTa$par
# 
# # Fit model for various samples
# ModelFitSamples_pracma <- lapply(1:10, function(x) {
#   
#   cat("\n********   Estimating effect of insertion length: sample", x, "   **********\n")
#   blnL1widthNA <- is.na(L1TotData$L1width)
#   L1TotData$L1width_sample <- NA
#   L1TotData$L1width_sample[!blnL1widthNA] <- SampleTrueL1Width(
#     SimL1widthTrue = L1DetectAgg_withL1$L1widthTrue, 
#     SimL1widthEst  = L1DetectAgg_withL1$L1widthEst,
#     SimL1Freq = L1DetectAgg_withL1$EstFreq,
#     EstL1width = L1TotData$L1width[!blnL1widthNA],
#     ObsL1Freq = L1TotData$Freq[!blnL1widthNA],
#     L1widthBreaks = seq(0, 6500, 500))
#    
#   # Create a matrix of predictor variables (L1 start and boolean variable for)
#   PredictMat <- L1TotData[, c("blnFull", "L1width_sample", "Freq")]
#   blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))
#   
#   FitSelectionModels_pracma(PredictMat[!blnNA, ],  
#                                   Freqs = round(L1TotData$Freq[!blnNA], 0), 
#                                   Counts = rep(1, sum(!blnNA)), 
#                                   PopSize = PopSize, 
#                                   SampleSize = L1TotData$SampleSize[!blnNA],
#                                   blnIns = L1TotData$blnIns[!blnNA], 
#                                   LogRegCoeff = LogRegL1Ref$coefficients,
#                                   DetectProb = L1TotData$DetectProb[!blnNA],
#                                   aBorder = 0.003, 
#                                   bBorder = 10^(-2), 
#                                   cBorder = 10^(-5))
#   
# })
# 
# # Summarize samples
# SampleSummary <- t(sapply(ModelFitSamples, function(x){
#   AIC       <- as.numeric(as.character(x$AICTab$AIC))
#   Predictor <- as.character(x$AICTab$Predictor)
#   Pars      <- as.character(x$AICTab$Pars)
#   SortedDiffs <- sort(AIC - min(AIC))
#   c(BestModelPredict =  Predictor[which.min(AIC)],
#     BestModelPars    =  Pars[which.min(AIC)],
#     DeltaAIC = SortedDiffs[2])
# }))

# Save everything
save.image(SelectResultOutPath)
load(SelectResultOutPath)

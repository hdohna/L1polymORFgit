# The script below estimates selection coefficients of L1 from the 
# 1000 genome data 
# 

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
DataPath            <- 'D:/L1polymORF/Data/'
MeltInsPath         <- "D:/L1polymORF/Data/nstd144.GRCh37.variant_call.vcf"
MeltDelPath         <- "D:/L1polymORF/Data/DEL.final_comp.vcf"
ChrLPath            <- 'D:/L1polymORF/Data/ChromLengthsHg19.Rdata'
InputPath           <- 'D:/L1polymORF/Data/SingletonAnalysis_unphased.RData'
L1RefPath           <- 'D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv'
L1RefRangePath      <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
RegrOutputPath      <- "D:/L1polymORF/Data/L1RegressionResults.RData"
SelectTabOutPath    <- "D:/L1polymORF/Data/L1SelectionResults_MELT.csv"
SelectGenTabOutPath <- "D:/L1polymORF/Data/L1SelectionGeneResults_MELT.csv"
SelectResultOutPath <- "D:/L1polymORF/Data/L1SelectionResults_MELT.RData"
SelectWithinGenTabOutPath <- "D:/L1polymORF/Data/L1SelectionWithinGeneResults_MELT.csv"
SelectSingletonTabOutPath <- "D:/L1polymORF/Data/L1SelectionSingletonResults_MELT.csv"

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

# Load simulation results
load("D:/L1polymORF/Data/L1Simulated_AdditionalInfo_MELT.RData")

# Read in vcf file with MELT insertion calls
MEInsCall <- read.table(MeltInsPath, 
                        as.is = T,
                        col.names = c("Chrom", "Pos", "ID", "Alt", "Type", "V6", 
                                      "V7", "Info"))
MEInsCall <- MEInsCall[MEInsCall$Type == "<INS:ME:LINE1>",]

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

# Load previously generated objects
load(InputPath)
load(L1GRPath)
load(ChrLPath)
load(L1RefRangePath)
load(RegrOutputPath)
load("D:/L1polymORF/Data/DelVsL1Length.RData")

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

# Number of L1 that are fixed at proportion 1
N1 <- length(L1GRanges) - length(OL_MEDelRefL1@from)

RefL1Data <- RefL1Data[OL_MEDelRefL1@from, ]
L1GRanges <- L1GRanges[OL_MEDelRefL1@from]

# Put data of non-reference L1 (insertions) and reference L1 (deletions) 
# together
L1TotData <- rbind(MEInsCall[ ,c("L1width", "Freq", "SampleSize", "blnFull")],
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

# Combine genomic ranges
L1TotGR <- c(MEIns_GR, L1GRanges)

# Number of L1 that are not fixed
Nnf <- nrow(L1TotData)

# Add sampled L1 width
L1DetectAgg_withL1 <- L1DetectAgg_Group[
  which(L1DetectAgg_Group$L1widthTrue > 0 & L1DetectAgg_Group$L1widthEst > 0 ),]
blnL1widthNA <- is.na(L1TotData$L1width)
L1TotData$L1width_sample1 <- NA
L1TotData$L1width_sample2 <- NA
L1TotData$L1width_sample3 <- NA
L1TotData$L1width_sample4 <- NA
L1TotData$L1width_sample1[!blnL1widthNA] <- SampleTrueL1Width(
    SimL1widthTrue = L1DetectAgg_withL1$L1widthTrue, 
    SimL1widthEst  = L1DetectAgg_withL1$L1widthEst,
    SimL1Freq = L1DetectAgg_withL1$EstFreq,
    EstL1width = L1TotData$L1width[!blnL1widthNA],
    ObsL1Freq = L1TotData$Freq[!blnL1widthNA],
    L1widthBreaks = seq(0, 6500, 500))
L1TotData$L1width_sample2[!blnL1widthNA] <- SampleTrueL1Width(
  SimL1widthTrue = L1DetectAgg_withL1$L1widthTrue, 
  SimL1widthEst  = L1DetectAgg_withL1$L1widthEst,
  SimL1Freq = L1DetectAgg_withL1$EstFreq,
  EstL1width = L1TotData$L1width[!blnL1widthNA],
  ObsL1Freq = L1TotData$Freq[!blnL1widthNA],
  L1widthBreaks = seq(0, 6500, 500))
L1TotData$L1width_sample3[!blnL1widthNA] <- SampleTrueL1Width(
  SimL1widthTrue = L1DetectAgg_withL1$L1widthTrue, 
  SimL1widthEst  = L1DetectAgg_withL1$L1widthEst,
  SimL1Freq = L1DetectAgg_withL1$EstFreq,
  EstL1width = L1TotData$L1width[!blnL1widthNA],
  ObsL1Freq = L1TotData$Freq[!blnL1widthNA],
  L1widthBreaks = seq(0, 6500, 500))
L1TotData$L1width_sample4[!blnL1widthNA] <- SampleTrueL1Width(
  SimL1widthTrue = L1DetectAgg_withL1$L1widthTrue, 
  SimL1widthEst  = L1DetectAgg_withL1$L1widthEst,
  SimL1Freq = L1DetectAgg_withL1$EstFreq,
  EstL1width = L1TotData$L1width[!blnL1widthNA],
  ObsL1Freq = L1TotData$Freq[!blnL1widthNA],
  L1widthBreaks = seq(0, 6500, 500))

cat("done!\n")

###################################################
#                                                 #
#   Fit effect of insertion length on selection   #
#                                                 #
###################################################

# Create a matrix of predictor variables (L1 start and boolean variable for)
PredictMat <- L1TotData[, c("blnFull", "L1width", "L1width_sample1", 
                            "L1width_sample2", "L1width_sample3", 
                            "L1width_sample4")]
blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))

cat("\n********   Estimating effect of insertion length: true length   **********\n")
ModelFit1 <- FitSelectionModels(PredictMat[!blnNA, 1:3],  
                                Freqs = round(L1TotData$Freq[!blnNA], 0), 
                                Counts = rep(1, sum(!blnNA)), 
                                PopSize = PopSize, 
                                SampleSize = L1TotData$SampleSize[!blnNA],
                                blnIns = L1TotData$blnIns[!blnNA], 
                                LogRegCoeff = LogRegL1Ref$coefficients,
                                DetectProb = L1TotData$DetectProb[!blnNA],
                                aBorder = 0.003, 
                                bBorder = 10^(-5), 
                                cBorder = 10^(-2))

cat("\n********   Estimating effect of insertion length: sampled length 1   **********\n")
ModelFit2 <- FitSelectionModels(PredictMat[!blnNA, c(1, 3, 4)],  
                                Freqs = round(L1TotData$Freq[!blnNA], 0), 
                                Counts = rep(1, sum(!blnNA)), 
                                PopSize = PopSize, 
                                SampleSize = L1TotData$SampleSize[!blnNA],
                                blnIns = L1TotData$blnIns[!blnNA], 
                                LogRegCoeff = LogRegL1Ref$coefficients,
                                DetectProb = L1TotData$DetectProb[!blnNA],
                                aBorder = 0.003, 
                                bBorder = 10^(-5), 
                                cBorder = 10^(-2))

cat("\n********   Estimating effect of insertion length: sampled length 2   **********\n")
ModelFit3 <- FitSelectionModels(PredictMat[!blnNA, c(1, 4, 5)],  
                                Freqs = round(L1TotData$Freq[!blnNA], 0), 
                                Counts = rep(1, sum(!blnNA)), 
                                PopSize = PopSize, 
                                SampleSize = L1TotData$SampleSize[!blnNA],
                                blnIns = L1TotData$blnIns[!blnNA], 
                                LogRegCoeff = LogRegL1Ref$coefficients,
                                DetectProb = L1TotData$DetectProb[!blnNA],
                                aBorder = 0.003, 
                                bBorder = 10^(-5), 
                                cBorder = 10^(-2))

cat("\n********   Estimating effect of insertion length: sampled length 3   **********\n")
ModelFit4 <- FitSelectionModels(PredictMat[!blnNA, c(1, 5, 6)],  
                                Freqs = round(L1TotData$Freq[!blnNA], 0), 
                                Counts = rep(1, sum(!blnNA)), 
                                PopSize = PopSize, 
                                SampleSize = L1TotData$SampleSize[!blnNA],
                                blnIns = L1TotData$blnIns[!blnNA], 
                                LogRegCoeff = LogRegL1Ref$coefficients,
                                DetectProb = L1TotData$DetectProb[!blnNA],
                                aBorder = 0.003, 
                                bBorder = 10^(-5), 
                                cBorder = 10^(-2))

cat("\n********   Estimating effect of insertion length: sampled length 4   **********\n")
ModelFit5 <- FitSelectionModels(PredictMat[!blnNA, c(1, 6, 2)],  
                                Freqs = round(L1TotData$Freq[!blnNA], 0), 
                                Counts = rep(1, sum(!blnNA)), 
                                PopSize = PopSize, 
                                SampleSize = L1TotData$SampleSize[!blnNA],
                                blnIns = L1TotData$blnIns[!blnNA], 
                                LogRegCoeff = LogRegL1Ref$coefficients,
                                DetectProb = L1TotData$DetectProb[!blnNA],
                                aBorder = 0.003, 
                                bBorder = 10^(-5), 
                                cBorder = 10^(-2))

# Save everything
save.image(SelectResultOutPath)

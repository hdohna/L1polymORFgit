# The script below estimates selection coefficients of L1 from the 
# 1000 genome data using insertion estimates obtained by MELT 
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
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

# Specify file paths
DataPath            <- 'D:/L1polymORF/Data/'
MeltInsPath         <- "D:/L1polymORF/Data/L1_SingleMELT_fullGenome_CombinedVcfs"
MeltDelPath         <- "D:/L1polymORF/Data/DEL.final_comp.vcf"
ChrLPath            <- 'D:/L1polymORF/Data/ChromLengthsHg19.Rdata'
InputPath           <- 'D:/L1polymORF/Data/SingletonAnalysis_unphased.RData'
L1RefPath           <- 'D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv'
L1RefRangePath      <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
RegrOutputPath      <- "D:/L1polymORF/Data/L1RegressionResults.RData"
SelectTabOutPath    <- "D:/L1polymORF/Data/L1SelectionResults_MELT_Single.csv"
SelectResultOutPath <- "D:/L1polymORF/Data/L1SelectionResults_MELT_Single.RData"
SampleInfoPath <- "D:/L1polymORF/Data/1000GenomeSampleInfo.txt"

# Specify logistic regression coefficients for relationship between insert
# size and detection probability
load("D:/L1polymORF/Data/L1Simulated_MELT.RData")
L1Detect$blnDetectPass <- L1Detect$blnDetect & L1Detect$EstFilter == "PASS"
LogReg_DetectL1width <- glm(blnDetectPass ~ L1widthTrue, 
                                data = L1Detect[L1Detect$L1GenoTrue == 1,],
                                family = binomial)
summary(LogReg_DetectL1width)
L1SizeDetectCoeff <- c(a = LogReg_DetectL1width$coefficients[1], 
                       b = LogReg_DetectL1width$coefficients[2])

# False discovery rate for selected L1
FDR <- 0.1

# Human effective population size
PopSize <- 10^5

# Minimum length for a full L1
MinLengthFullL1 <- 6000

# Sample size for ME insertion calls

# Coefficients for the probability to be reference, depending on L1 frequency
LogRegL1RefCoeff <- c(-4.706573, 9.737618)

##########################################
#                                        #
#     Define functions                   #
#                                        #
##########################################

# Function to create MEInsCall data
# CreateMEInsCall <- function(Samples2use, MEInsCallPerL1 = MEInsCallPerL1){
#   
#   MEInsCallPerL1 <- MEInsCallPerL1[MEInsCallPerL1$SampleID %in% Samples2use, ]
#   MEInsCall <- aggregate(MEInsCallPerL1$GenoNum, 
#                          by = list(MEInsCallPerL1$L1ID), 
#                          FUN = sum)
#   colnames(MEInsCall) <- c("L1ID", "Freq")
#   MEInsCall$CHROM <- sapply(MEInsCall$L1ID, function(x) strsplit(x, " ")[[1]][1])
#   MEInsCall$POS   <- sapply(MEInsCall$L1ID, function(x) as.numeric(strsplit(x, " ")[[1]][2]))
#   L1IDmatch       <- match(MEInsCall$L1ID, MEInsCallPerL1$L1ID)
#   MEInsCall$L1width <- MEInsCallPerL1$L1width[L1IDmatch]
#   
#   # Add columns necessary for analysis 
#   MEInsCall$AF <- MEInsCall$Freq / MEInsSamplesize
#   MEInsCall$SampleSize <- 2 * MEInsSamplesize
#   MEInsCall$blnFull    <- MEInsCall$L1width >= MinLengthFullL1
#   MEInsCall
# }


##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("\n\nLoading and processing data ...\n")

SampleInfo <- read.table(SampleInfoPath, header = T)

# Read in vcf file with MELT insertion calls
MEInsCallPerL1  <- read.table(MeltInsPath, as.is = T)
NrSamples       <- length(unique(MEInsCallPerL1$SampleID))
MEInsSamplesize <- 2*NrSamples

# Subset to get entries with genotype
idxWithGeno    <- which(nchar(MEInsCallPerL1$Genotype) > 0)
MEInsCallPerL1 <- MEInsCallPerL1[idxWithGeno, ]

# Add columns for numeric genotype, insertion length
cat("Adding L1 info ...")
MEInsCallPerL1$GenoNum <- sapply(MEInsCallPerL1$Genotype, GetFromVcfGeno_GenoNum)
MEInsCallPerL1$L1width <- sapply(MEInsCallPerL1$INFO, GetFromVcfINFO_SVLength)
MEInsCallPerL1$PropCovered  <- sapply(MEInsCallPerL1$INFO, 
                                      GetFromVcfINFO_MELT_PropCovered)
MEInsCallPerL1$L1Diff  <- sapply(MEInsCallPerL1$INFO, GetFromVcfINFO_MELT_PropDiff)
L1StartEnd  <- t(sapply(MEInsCallPerL1$INFO, GetFromVcfINFO_MELT_L1StartEnd))
colnames(L1StartEnd) <- c("L1Start", "L1End")
MEInsCallPerL1 <- cbind(MEInsCallPerL1, L1StartEnd)
cat("done!\n")

# Get L1 ID and a
MEInsCallPerL1$L1ID <- CollapseClosePos_idx(MEInsCallPerL1, 
                                            ChromCol = "X.CHROM", 
                                            PosCol = "POS",
                                            OLRange = 10)

# Aggregate values per L1  
cat("Aggregating data per L1 ...")
MEInsCall <- AggDataFrame(MEInsCallPerL1, GroupCol = "L1ID", 
                          MeanCols = "L1Diff", 
                          SumCols = "GenoNum",
                          MedCols = "L1width", 
                          LengthCols = "L1Diff",
                          MinCols = "L1Start",
                          MaxCols = "L1End", Addcols = c("X.CHROM", "POS"))
colnames(MEInsCall)[colnames(MEInsCall) == "GenoNum_sum"] <- "Freq"
colnames(MEInsCall)[colnames(MEInsCall) == "X.CHROM"]     <- "CHROM"
colnames(MEInsCall)[colnames(MEInsCall) == "L1width_med"] <- "L1width"

cat("done!\n")

# Add columns necessary for analysis 
MEInsCall$AF         <- MEInsCall$Freq / MEInsSamplesize
MEInsCall$SampleSize <- 2 * MEInsSamplesize
MEInsCall$blnFull    <- MEInsCall$L1width >= MinLengthFullL1

# Create GRanges object for MEInsCall
MEInsCall$ChromName <- paste("chr", MEInsCall$CHROM, sep = "")
MEIns_GR <- makeGRangesFromDataFrame(df = MEInsCall,
                                     seqnames.field = "ChromName",
                                     start.field = "Pos",
                                     end.field = "Pos")



# Read in vcf file with MELT insertion calls
MEInsCall_MELT <- read.table("D:/L1polymORF/Data/nstd144.GRCh37.variant_call.vcf", 
                        as.is = T,
                        col.names = c("Chrom", "Pos", "ID", "Alt", "Type", "V6", 
                                      "V7", "Info"))
MEInsCall_MELT <- MEInsCall_MELT[MEInsCall_MELT$Type == "<INS:ME:LINE1>",]
MEInsCall_MELT$ChromName <- paste("chr", MEInsCall_MELT$Chrom, sep = "")
MEIns_GR_MELT <- makeGRangesFromDataFrame(df = MEInsCall_MELT,
                                     seqnames.field = "ChromName",
                                     start.field = "Pos",
                                     end.field = "Pos")
DistObj <- Dist2Closest(MEIns_GR, MEIns_GR_MELT)
DistObj

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
L1TotData$DetectProb[L1TotData$blnIns] <- 
  exp(L1SizeDetectCoeff[1] + 
        L1SizeDetectCoeff[2] * L1TotData$L1width[L1TotData$blnIns]) / 
  (1 + exp(L1SizeDetectCoeff[1] + 
             L1SizeDetectCoeff[2]*L1TotData$L1width[L1TotData$blnIns]))

# Remove all L1 with zero frequency
L1TotData[L1TotData$Freq == 0, ]
L1TotData <- L1TotData[L1TotData$Freq > 0, ]
hist(L1TotData$Freq, breaks = 0:1000, xlim = c(0, 10))

# Perform logistic regression for the probability of reference L1 as function
# of L1 frequency
L1TotData$blnRef <- !L1TotData$blnIns
LogRegL1Ref <- glm(blnRef ~ L1Freq, family = binomial, data = L1TotData)
LogRegL1Ref$coefficients

# Combine genomic ranges
L1TotGR <- c(MEIns_GR, L1GRanges)


# Number of L1 that are not fixed
Nnf <- nrow(L1TotData)

# Read information about 1000 genome samples
SampleInfo  <- read.table(G1000SamplePath, header = T)
SampleMatch <- match(SampleColumns, SampleInfo$sample)
Pops        <- SampleInfo$super_pop[SampleMatch]
NrS         <- length(SampleColumns)

cat("done!\n")


###################################################
#                                                 #
#   Fit effect of insertion length on selection   #
#                                                 #
###################################################

cat("\n********   Estimating effect of insertion length    **********\n")

# Create a matrix of predictor variables (L1 start and boolean variable for)
PredictMat <- L1TotData[, c("L1width", "blnFull", "Freq", "SampleSize", "blnIns")]
                        
blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,]))) |
  L1TotData$Freq == 0
sum(!blnNA)
which(L1TotData$Freq == 0)
max(L1TotData$Freq / L1TotData$SampleSize, na.rm = T)
max(L1TotData$Freq, na.rm = T)

# Estimate maximum likelihood for a single selection coefficient
cat("Estimate maximum likelihood for a single selection coefficient\n")
ML_1Par <-  constrOptim(theta = c(a = 0),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = round(L1TotData$Freq[!blnNA], 0),
                            Counts = rep(1, sum(!blnNA)),
                            Predict = PredictMat[!blnNA, 1:3],
                            a = x[1], b = 0, c = 0, d = 0, N = PopSize,
                            SampleSize = L1TotData$SampleSize[!blnNA],
                            blnIns = L1TotData$blnIns[!blnNA], 
                            LogRegCoeff = LogRegL1Ref$coefficients,
                            DetectProb = L1TotData$DetectProb[!blnNA]),
                          grad = NULL,
                          ui = rbind(1,-1),
                          ci = c(a = -0.003, a = -0.003),
                          method = "Nelder-Mead")
cat("done!\n")


# Get maximum likelihood estimate for effect of L1 start on selection
cat("Estimate effect of L1 start on selections ...")
ML_L1width <-  constrOptim(theta = c(a = ML_1Par$par, b = 0),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = round(L1TotData$Freq[!blnNA], 0),
                            Counts = rep(1, sum(!blnNA)),
                            Predict = PredictMat[!blnNA, 1:3],
                            a = x[1], b = x[2], c = 0, d = 0, N = PopSize, 
                            SampleSize = L1TotData$SampleSize[!blnNA],
                            blnIns = L1TotData$blnIns[!blnNA], 
                            LogRegCoeff = LogRegL1Ref$coefficients,
                            DetectProb = L1TotData$DetectProb[!blnNA]),
                          grad = NULL,
                          ui = rbind(c(1, 0),  c(0, 1),   
                                     c(-1, 0), c(0, -1)),
                          ci = c(a = -0.02, c = -10^(-6), 
                                 a = -0.02, c = -10^(-6)),
                          method = "Nelder-Mead")
cat("done!\n")

# Get maximum likelihood estimate for effect of full-length L1 on selection
cat("Estimate effect of L1 full-length on selections ...")
ML_L1full <-  constrOptim(theta = c(a = ML_1Par$par, c = 0),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = round(L1TotData$Freq[!blnNA], 0),
                            Counts = rep(1, sum(!blnNA)),
                            Predict = PredictMat[!blnNA, 1:3],
                            a = x[1], b = 0, c = x[2], d = 0, N = PopSize, 
                            SampleSize = L1TotData$SampleSize[!blnNA],
                            blnIns = L1TotData$blnIns[!blnNA], 
                            LogRegCoeff = LogRegL1Ref$coefficients,
                            DetectProb = L1TotData$DetectProb[!blnNA]),
                          grad = NULL,
                      ui = rbind(c(1, 0),  c(0, 1),   
                                 c(-1, 0), c(0, -1)),
                      ci = c(a = -0.02, d = -10^(-3), 
                             a = -0.02, d = -10^(-3)),
                      method = "Nelder-Mead")
cat("done!\n")

# Determine maximum likelihood with 3 parameters (selection coefficient as 
# function of L1 start and indicator for full-length)
cat("Maximizing likelihood for three parameters ...")
ML_L1widthL1full <- constrOptim(theta = c(a = ML_L1width$par[1], 
                                          b = ML_L1width$par[2], 
                                          c = ML_L1full$par[2]),
                  f = function(x) -AlleleFreqLogLik_4Par(
                    Freqs = round(L1TotData$Freq[!blnNA], 0),
                    Counts = rep(1, sum(!blnNA)),
                    Predict = PredictMat[!blnNA, 1:3],
                    a = x[1], b = x[2], c = x[3], d = 0, N = PopSize, 
                    SampleSize = L1TotData$SampleSize[!blnNA],
                    blnIns = L1TotData$blnIns[!blnNA], 
                    LogRegCoeff = LogRegL1Ref$coefficients,
                    DetectProb = L1TotData$DetectProb[!blnNA]),
                  grad = NULL,
                  ui = rbind(c(1, 0, 0),  c(0, 1, 0),  c(0, 0, 1), 
                             c(-1, 0, 0), c(0, -1, 0), c(0, 0, -1)),
                  ci = c(a = -0.01, b = -10^(-6), d = -10^(-3), 
                         a = -0.02, b = -10^(-6), d = -10^(-3)),
                  method = "Nelder-Mead")
cat("done!\n")

###################################################
#                                                 #
#  Summarize results                              #
#                                                 #
###################################################

# Function to extract AIC from optim results
GetAIC <- function(OptimResults){
  round(2 * (length(OptimResults$par) + OptimResults$value), 2)
}
GetParVals <- function(OptimResults){
  Results <- paste(names(OptimResults$par), 
                   format(OptimResults$par, digits = 2), sep = " = ",
                   collapse = ", ")
}
GetNPar <- function(OptimResults){
  length(OptimResults$par)
}

# Get columns of AIC and parameter values
Cols2Append <- t(sapply(list(ML_1Par, 
                             ML_L1width, 
                             ML_L1full, 
                             ML_L1widthL1full), function(x){
           c(AIC = GetAIC(x), Pars = GetParVals(x))
         }))
# Combine AIC values into one vector
AICTab <- cbind(data.frame(
            NrParameters = c(1, 2, 2, 3),
            Predictor = c("none", 
                          "L1 width", 
                          "L1 full-length", 
                          "L1 width and full-length"),
            stringsAsFactors = F),
            Cols2Append)
                     
# Save table with AIC
write.csv(AICTab, SelectTabOutPath)
save.image(SelectResultOutPath)

nrow(PredictMat)
ModelFit1 <- FitSelectionModels(PredictMat[!blnNA, 1:3],  
                               Freqs = round(L1TotData$Freq[!blnNA], 0), 
                               Counts = rep(1, sum(!blnNA)), 
                               PopSize = PopSize, 
                               SampleSize = L1TotData$SampleSize[!blnNA],
                               blnIns = L1TotData$blnIns[!blnNA], 
                               LogRegCoeff = LogRegL1Ref$coefficients,
                               DetectProb = L1TotData$DetectProb[!blnNA],
                               aBorder = 0.003, 
                               bBorder = 10^(-6), 
                               cBorder = 10^(-6)
                               
)

###################################################
#                                                 #
#  Same as above                             #
#                                                 #
###################################################



###################################################
#                                                        #
#   Fit effect of insertion length on selection,
#   only using full-length and fragments below 500 bp    #
#                                                        #
###################################################

# cat("\n********   Estimating effect of insertion length    **********\n")
# 
# # Create a matrix of predictor variables (L1 start and boolean variable for)
# blnSubset  <- (L1TotData$L1width <= 500 | L1TotData$L1width >= 6000) & (!is.na(L1TotData$Freq))
# PredictMat <- L1TotData[, c("L1width", "blnFull", "Freq", 
#                                             "SampleSize", "blnIns")]
# 
# blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,]))) |
#   L1TotData$Freq == 0 | (!blnSubset)
# sum(blnNA)
# which(L1TotData$Freq == 0)
# max(L1TotData$Freq / L1TotData$SampleSize, na.rm = T)
# max(L1TotData$Freq, na.rm = T)
# 
# # Estimate maximum likelihood for a single selection coefficient
# cat("Estimate maximum likelihood for a single selection coefficient\n")
# ML_1Par <-  constrOptim(theta = c(a = 0),
#                         f = function(x) -AlleleFreqLogLik_4Par(
#                           Freqs = round(L1TotData$Freq[!blnNA], 0),
#                           Counts = rep(1, sum(!blnNA)),
#                           Predict = PredictMat[!blnNA, 1:3],
#                           a = x[1], b = 0, c = 0, d = 0, N = PopSize,
#                           SampleSize = L1TotData$SampleSize[!blnNA],
#                           blnIns = L1TotData$blnIns[!blnNA], 
#                           LogRegCoeff = LogRegL1RefCoeff,
#                           DetectProb = L1TotData$DetectProb[!blnNA]),
#                         grad = NULL,
#                         ui = rbind(1,-1),
#                         ci = c(a = -0.003, a = -0.003),
#                         method = "Nelder-Mead")
# cat("done!\n")
# 
# 
# # Get maximum likelihood estimate for effect of L1 start on selection
# cat("Estimate effect of L1 start on selections ...")
# ML_L1width <-  constrOptim(theta = c(a = ML_1Par$par, b = 0),
#                            f = function(x) -AlleleFreqLogLik_4Par(
#                              Freqs = round(L1TotData$Freq[!blnNA], 0),
#                              Counts = rep(1, sum(!blnNA)),
#                              Predict = PredictMat[!blnNA, 1:3],
#                              a = x[1], b = x[2], c = 0, d = 0, N = PopSize, 
#                              SampleSize = L1TotData$SampleSize[!blnNA],
#                              blnIns = L1TotData$blnIns[!blnNA], 
#                              LogRegCoeff = LogRegL1RefCoeff,
#                              DetectProb = L1TotData$DetectProb[!blnNA]),
#                            grad = NULL,
#                            ui = rbind(c(1, 0),  c(0, 1),   
#                                       c(-1, 0), c(0, -1)),
#                            ci = c(a = -0.02, c = -10^(-6), 
#                                   a = -0.02, c = -10^(-6)),
#                            method = "Nelder-Mead")
# cat("done!\n")
# 
# # Get maximum likelihood estimate for effect of full-length L1 on selection
# cat("Estimate effect of L1 full-length on selections ...")
# ML_L1full <-  constrOptim(theta = c(a = ML_1Par$par, c = 0),
#                           f = function(x) -AlleleFreqLogLik_4Par(
#                             Freqs = round(L1TotData$Freq[!blnNA], 0),
#                             Counts = rep(1, sum(!blnNA)),
#                             Predict = PredictMat[!blnNA, 1:3],
#                             a = x[1], b = 0, c = x[2], d = 0, N = PopSize, 
#                             SampleSize = L1TotData$SampleSize[!blnNA],
#                             blnIns = L1TotData$blnIns[!blnNA], 
#                             LogRegCoeff = LogRegL1RefCoeff,
#                             DetectProb = L1TotData$DetectProb[!blnNA]),
#                           grad = NULL,
#                           ui = rbind(c(1, 0),  c(0, 1),   
#                                      c(-1, 0), c(0, -1)),
#                           ci = c(a = -0.02, d = -10^(-3), 
#                                  a = -0.02, d = -10^(-3)),
#                           method = "Nelder-Mead")
# cat("done!\n")
# 
# # Determine maximum likelihood with 3 parameters (selection coefficient as 
# # function of L1 start and indicator for full-length)
# cat("Maximizing likelihood for three parameters ...")
# ML_L1widthL1full <- constrOptim(theta = c(a = ML_L1width$par[1], 
#                                           b = ML_L1width$par[2], 
#                                           c = ML_L1full$par[2]),
#                                 f = function(x) -AlleleFreqLogLik_4Par(
#                                   Freqs = round(L1TotData$Freq[!blnNA], 0),
#                                   Counts = rep(1, sum(!blnNA)),
#                                   Predict = PredictMat[!blnNA, 1:3],
#                                   a = x[1], b = x[2], c = x[3], d = 0, N = PopSize, 
#                                   SampleSize = L1TotData$SampleSize[!blnNA],
#                                   blnIns = L1TotData$blnIns[!blnNA], 
#                                   LogRegCoeff = LogRegL1RefCoeff,
#                                   DetectProb = L1TotData$DetectProb[!blnNA]),
#                                 grad = NULL,
#                                 ui = rbind(c(1, 0, 0),  c(0, 1, 0),  c(0, 0, 1), 
#                                            c(-1, 0, 0), c(0, -1, 0), c(0, 0, -1)),
#                                 ci = c(a = -0.01, b = -10^(-6), d = -10^(-3), 
#                                        a = -0.02, b = -10^(-6), d = -10^(-3)),
#                                 method = "Nelder-Mead")
# cat("done!\n")
# 
# ###################################################
# #                                                 #
# #  Summarize results                              #
# #                                                 #
# ###################################################
# 
# # Function to extract AIC from optim results
# GetAIC <- function(OptimResults){
#   round(2 * (length(OptimResults$par) + OptimResults$value), 2)
# }
# GetParVals <- function(OptimResults){
#   Results <- paste(names(OptimResults$par), 
#                    format(OptimResults$par, digits = 2), sep = " = ",
#                    collapse = ", ")
# }
# GetNPar <- function(OptimResults){
#   length(OptimResults$par)
# }
# 
# # Get columns of AIC and parameter values
# Cols2Append <- t(sapply(list(ML_1Par, 
#                              ML_L1width, 
#                              ML_L1full, 
#                              ML_L1widthL1full), function(x){
#                                c(AIC = GetAIC(x), Pars = GetParVals(x))
#                              }))
# # Combine AIC values into one vector
# AICTab <- cbind(data.frame(
#   NrParameters = c(1, 2, 2, 3),
#   Predictor = c("none", 
#                 "L1 width", 
#                 "L1 full-length", 
#                 "L1 width and full-length"),
#   stringsAsFactors = F),
#   Cols2Append)
# 
# # Save table with AIC
# write.csv(AICTab, SelectTabOutPath)
# save.image(SelectResultOutPath)
# 
###################################################
#                                                 #
#   Fit effect of insertion length on selection
#   without length-dependent detection            #
#                                                 #
###################################################

cat("\n********   Estimating effect of insertion length    **********\n")

# # Create a matrix of predictor variables (L1 start and boolean variable for)
# PredictMat <- L1TotData[, c("L1width", "blnFull", "Freq", "SampleSize", "blnIns")]
# 
# blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,]))) |
#   L1TotData$Freq == 0
# sum(!blnNA)
# which(L1TotData$Freq == 0)
# max(L1TotData$Freq / L1TotData$SampleSize, na.rm = T)
# max(L1TotData$Freq, na.rm = T)
# 
# # Estimate maximum likelihood for a single selection coefficient
# cat("Estimate maximum likelihood for a single selection coefficient\n")
# ML_1Par <-  constrOptim(theta = c(a = 0),
#                         f = function(x) -AlleleFreqLogLik_4Par(
#                           Freqs = round(L1TotData$Freq[!blnNA], 0),
#                           Counts = rep(1, sum(!blnNA)),
#                           Predict = PredictMat[!blnNA, 1:3],
#                           a = x[1], b = 0, c = 0, d = 0, N = PopSize,
#                           SampleSize = L1TotData$SampleSize[!blnNA],
#                           blnIns = L1TotData$blnIns[!blnNA], 
#                           LogRegCoeff = LogRegL1RefCoeff,
#                           DetectProb = rep(mean(L1TotData$DetectProb, na.rm = T),
#                                            sum(!blnNA))),
#                         grad = NULL,
#                         ui = rbind(1,-1),
#                         ci = c(a = -0.003, a = -0.003),
#                         method = "Nelder-Mead")
# cat("done!\n")
# 
# 
# # Get maximum likelihood estimate for effect of L1 start on selection
# cat("Estimate effect of L1 start on selections ...")
# ML_L1width <-  constrOptim(theta = c(a = ML_1Par$par, b = 0),
#                            f = function(x) -AlleleFreqLogLik_4Par(
#                              Freqs = round(L1TotData$Freq[!blnNA], 0),
#                              Counts = rep(1, sum(!blnNA)),
#                              Predict = PredictMat[!blnNA, 1:3],
#                              a = x[1], b = x[2], c = 0, d = 0, N = PopSize, 
#                              SampleSize = L1TotData$SampleSize[!blnNA],
#                              blnIns = L1TotData$blnIns[!blnNA], 
#                              LogRegCoeff = LogRegL1RefCoeff,
#                              DetectProb = rep(mean(L1TotData$DetectProb, na.rm = T),
#                                               sum(!blnNA))),                           
#                            grad = NULL,
#                            ui = rbind(c(1, 0),  c(0, 1),   
#                                       c(-1, 0), c(0, -1)),
#                            ci = c(a = -0.02, c = -10^(-6), 
#                                   a = -0.02, c = -10^(-6)),
#                            method = "Nelder-Mead")
# cat("done!\n")
# 
# # Get maximum likelihood estimate for effect of full-length L1 on selection
# cat("Estimate effect of L1 full-length on selections ...")
# ML_L1full <-  constrOptim(theta = c(a = ML_1Par$par, c = 0),
#                           f = function(x) -AlleleFreqLogLik_4Par(
#                             Freqs = round(L1TotData$Freq[!blnNA], 0),
#                             Counts = rep(1, sum(!blnNA)),
#                             Predict = PredictMat[!blnNA, 1:3],
#                             a = x[1], b = 0, c = x[2], d = 0, N = PopSize, 
#                             SampleSize = L1TotData$SampleSize[!blnNA],
#                             blnIns = L1TotData$blnIns[!blnNA], 
#                             LogRegCoeff = LogRegL1RefCoeff,
#                             DetectProb = rep(mean(L1TotData$DetectProb, na.rm = T),
#                                              sum(!blnNA))),                          
#                           grad = NULL,
#                           ui = rbind(c(1, 0),  c(0, 1),   
#                                      c(-1, 0), c(0, -1)),
#                           ci = c(a = -0.02, d = -10^(-3), 
#                                  a = -0.02, d = -10^(-3)),
#                           method = "Nelder-Mead")
# cat("done!\n")
# 
# # Determine maximum likelihood with 3 parameters (selection coefficient as 
# # function of L1 start and indicator for full-length)
# cat("Maximizing likelihood for three parameters ...")
# ML_L1widthL1full <- constrOptim(theta = c(a = ML_L1width$par[1], 
#                                           b = ML_L1width$par[2], 
#                                           c = ML_L1full$par[2]),
#                                 f = function(x) -AlleleFreqLogLik_4Par(
#                                   Freqs = round(L1TotData$Freq[!blnNA], 0),
#                                   Counts = rep(1, sum(!blnNA)),
#                                   Predict = PredictMat[!blnNA, 1:3],
#                                   a = x[1], b = x[2], c = x[3], d = 0, N = PopSize, 
#                                   SampleSize = L1TotData$SampleSize[!blnNA],
#                                   blnIns = L1TotData$blnIns[!blnNA], 
#                                   LogRegCoeff = LogRegL1RefCoeff,
#                                   DetectProb = rep(mean(L1TotData$DetectProb, na.rm = T),
#                                                    sum(!blnNA))),                                
#                                 grad = NULL,
#                                 ui = rbind(c(1, 0, 0),  c(0, 1, 0),  c(0, 0, 1), 
#                                            c(-1, 0, 0), c(0, -1, 0), c(0, 0, -1)),
#                                 ci = c(a = -0.01, b = -10^(-6), d = -10^(-3), 
#                                        a = -0.02, b = -10^(-6), d = -10^(-3)),
#                                 method = "Nelder-Mead")
# cat("done!\n")
# 
# ###################################################
# #                                                 #
# #  Summarize results                              #
# #                                                 #
# ###################################################
# 
# 
# # Get columns of AIC and parameter values
# Cols2Append <- t(sapply(list(ML_1Par, 
#                              ML_L1width, 
#                              ML_L1full, 
#                              ML_L1widthL1full), function(x){
#                                c(AIC = GetAIC(x), Pars = GetParVals(x))
#                              }))
# # Combine AIC values into one vector
# AICTab <- cbind(data.frame(
#   NrParameters = c(1, 2, 2, 3),
#   Predictor = c("none", 
#                 "L1 width", 
#                 "L1 full-length", 
#                 "L1 width and full-length"),
#   stringsAsFactors = F),
#   Cols2Append)
# 
# # Save table with AIC
# write.csv(AICTab, SelectTabOutPath)
# save.image(SelectResultOutPath)

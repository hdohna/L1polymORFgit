# The script below estimates selection coefficients of L1 from the 
# 1000 genome data accounting for variable L1 length estimates 

##########################################
#                                        #
#     Load packages                      #
#                                        #
##########################################

# Source start script
source('D:/OneDrive - American University of Beirut/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(GenomicRanges)
library(pracma)

##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

# Specify file paths
DataPath            <- 'D:/OneDrive - American University of Beirut/L1polymORF/Data/'
MeltInsPath         <- "D:/OneDrive - American University of Beirut/L1polymORF/Data/nstd144.GRCh37.variant_call.vcf"
MeltDelPath         <- "D:/OneDrive - American University of Beirut/L1polymORF/Data/DEL.final_comp.vcf"
ChrLPath            <- 'D:/OneDrive - American University of Beirut/L1polymORF/Data/ChromLengthsHg19.Rdata'
InputPath           <- 'D:/OneDrive - American University of Beirut/L1polymORF/Data/SingletonAnalysis_unphased.RData'
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

# Breaks for L1 length classes
L1widthBreaks <- seq(0, 6500, 500)
L1MidPts <- 0.5*(L1widthBreaks[-1] + L1widthBreaks[-length(L1widthBreaks)])

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("\n\nLoading and processing data ...")

# Load simulation results
load("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1Simulated_AdditionalInfo_MELT.RData")

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
L1GRPath <- gsub("D:/", "D:/OneDrive - American University of Beirut/", L1GRPath)
load(L1GRPath)
ChrLPath <- gsub("D:/", "D:/OneDrive - American University of Beirut/", ChrLPath)
load(ChrLPath)
L1RefRangePath <- gsub("D:/", "D:/OneDrive - American University of Beirut/", 
                       L1RefRangePath)
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

# Group L1 width into classes
L1TotData$L1widthClass <- cut(L1TotData$L1width, breaks = L1widthBreaks)

# Group L1 width into classes
# L1WidthFreqTable <- table(L1TotData$L1widthClass, L1TotData$Freq)
# L1WidthFreqDF1 <- as.data.frame(L1WidthFreqTable)
# L1WidthFreqDF1 <- L1WidthFreqDF1[L1WidthFreqDF1$Var2 != "0", ]
L1WidthFreqDF <- aggregate(rep(1, nrow(L1TotData)), 
                           by = list(L1TotData$L1widthClass, L1TotData$Freq),
              FUN = sum)
L1WidthFreqDF <- L1WidthFreqDF[L1WidthFreqDF$Group.2 > 0, ]
sum(table(L1DetectAgg_Group$L1widthTrue, L1DetectAgg_Group$L1widthEst) > 0)
length(unique(L1DetectAgg_Group$L1widthEst)) * length(unique(L1DetectAgg_Group$L1widthEst))

# Subset to have only entries with non
L1DetectAgg_withL1 <- L1DetectAgg_Group[
  which(L1DetectAgg_Group$L1widthTrue > 0 & L1DetectAgg_Group$L1widthEst > 0 ),]
L1DetectAgg_withL1$L1widthTrue
L1DetectAgg_withL1

# L1SampleMat <- SampleMatTrueL1Width(SimL1widthTrue = L1DetectAgg_withL1$L1widthTrue, 
#                        SimL1widthEst  = L1DetectAgg_withL1$L1widthEst,
#                        SimL1Freq = L1DetectAgg_withL1$EstFreq,
#                        EstL1width = L1TotData$L1width[!blnL1widthNA],
#                        ObsL1Freq = L1TotData$Freq[!blnL1widthNA],
#                        L1widthBreaks = L1widthBreaks)

# Cut true and estimated L1 width into bins
SimL1widthTrueCut <- cut(L1DetectAgg_withL1$L1widthTrue, breaks = L1widthBreaks)
SimL1widthEstCut  <- cut(L1DetectAgg_withL1$L1widthEst, breaks  = L1widthBreaks)

# Create boolean indicator for estimated and true L1 width to be in the same 
# bin
blnSameBin <- SimL1widthTrueCut == SimL1widthEstCut

# Fit logistic regression for probability of estimated length in same bin as
# true length against estimated frequency
LogReg <- glm(blnSameBin ~ SimL1Freq, family = binomial)
summary(LogReg)
a <- LogReg$coefficients[1]
b <- LogReg$coefficients[2]
ExpPart  <- exp(LogReg$coefficients[1] + LogReg$coefficients[2] * ObsL1Freq)
ProbSame <- ExpPart / (1 + ExpPart)
plot(ObsL1Freq, ProbSame, col = rgb(0, 0, 0, 0.1), pch = 16)

# Get combinations of L1 width classes among estimated 
L1combos <- table(SimL1widthEstCut, SimL1widthTrueCut)
diag(L1combos) <- 0

# Matrix of true L1 width given estimated width
L1TrueGivenEst <- L1combos / rowSums(L1combos)

# List that gives for each estimated L1 length the probabilities of true
# L1 lengths and the lengths
L1TrueGivenEstList <- lapply(1:nrow(L1TrueGivenEst), function(x){
  idxNonZero <- which(L1TrueGivenEst[x, ] > 0)
  data.frame(idxNonZero = idxNonZero,
             Probs = L1TrueGivenEst[x, idxNonZero],
             L1Length = L1MidPts[idxNonZero])
})

# Add probability to be in the same length class
a <- LogReg$coefficients[1]
b <- LogReg$coefficients[2]
ExpPart  <- exp(LogReg$coefficients[1] + LogReg$coefficients[2] * L1TotData$Freq)
L1TotData$ProbSameLength <- ExpPart / (1 + ExpPart)

# Create dataframe with one row per frequency-length class combination
UniqueLengthClasses <- sort(unique(L1TotData$L1widthClass))
UniqueFreqs         <- sort(unique(L1TotData$Freq))
L1TotData$L1widthIdx <- match(L1TotData$L1widthClass, UniqueLengthClasses)
L1TotData$FreqIdx    <- match(L1TotData$Freq, UniqueFreqs)
L1TotData$FreqLength <- paste(L1TotData$Freq, L1TotData$L1widthClass, sep = "-")
L1TotData$N  <- 1
L1TotData$N  <- 1
FreqLengthDF <- AggDataFrame (L1TotData[L1TotData$blnIns,], 
                              GroupCol = "FreqLength", 
                              MeanCols = c("Freq", "L1widthIdx", "FreqIdx",
                                           "ProbSameLength"), 
                                         SumCols = "N")
blnNA <- sapply(1:nrow(FreqLengthDF), function(x) any(is.na(FreqLengthDF[x,])))
FreqLengthDF <- FreqLengthDF[!blnNA, ]

# Specify parameters that will go into function
a = 0
b = 0
c = 0
N = 10^5
LogRegCoeff = LogRegL1Ref$coefficients
DetectProb = 0.8
LengthFull = 6000
SampleSize = 2*MEInsSamplesize

# Determine which length classes are full-length
blnClassFull <- L1MidPts >= LengthFull

# Calculate selection coefficients
sVals = a + b*L1MidPts + c * blnClassFull

# Calculate integration constants
IntConsts <- sapply(sVals, function(x){
  AlleleFreqIntConst(s = x, N = N, SampleSize = SampleSize,
                     DetectProb = DetectProb,
                     LogRegCoeff = LogRegCoeff, blnIns = T)
})


i <- 1

LProbs <- sapply(1:nrow(FreqLengthDF), function(i){
  # Get index of current length
  idxLength <- FreqLengthDF$L1widthIdx_mean[i]
  
  # Get length values associated with current length and their characteristics
  LengthTrueEst <- L1TrueGivenEstList[[idxLength]] 
  blnFullOther <- LengthTrueEst$L1Length >= LengthFull
  sValsOther    <- a + b * LengthTrueEst$L1Length + c * blnFullOther
  
  # Calculate the probability of observed length at oberved frequency
  ProbObs <- FreqLengthDF$ProbSameLength_mean[i]  * 
    AlleleFreqSampleVarDiscrete1(FreqLengthDF$Freq_mean[i], 
                                 s = sVals[idxLength], 
                                 N = N, IntConst = IntConsts[idxLength], 
                                 SampleSize = SampleSize, DetectProb = DetectProb,
                                 LogRegCoeff = LogRegCoeff, blnIns = T) +
    
    (1 - FreqLengthDF$ProbSameLength_mean[i])  * 
    sum(sapply(1:nrow(LengthTrueEst), function(x){
      LengthTrueEst$Probs[x] *
        AlleleFreqSampleVarDiscrete1(FreqLengthDF$Freq_mean[i], 
                                     s = sValsOther[x], 
                                     N = N, IntConst = IntConsts[LengthTrueEst$idxNonZero[x]], 
                                     SampleSize = MEInsSamplesize, DetectProb = DetectProb,
                                     LogRegCoeff = LogRegCoeff, blnIns = T)
    }))
  
  FreqLengthDF$N_sum[i] * log(ProbObs)
  
})
sum(LProbs)

# Estimate maximum likelihood for a single selection coefficient
cat("Estimate maximum likelihood for a single selection coefficient\n")
ML_a <-  constrOptim(theta = c(a = 0),
                     f = function(x) 
                       -AlleleFreqLogLikVar_4Par(FreqLengthDF = FreqLengthDF, 
                                                 L1TrueGivenEstList = L1TrueGivenEstList, 
                                                 L1MidPts = L1MidPts,
                                                 a = x[1], 
                                                 b = 0, 
                                                 c = 0, 
                                                 d = 0, 
                                                 SD = NULL, 
                                                 N = PopSize, 
                                                 SampleSize = SampleSize,
                                                 blnIns = T, 
                                                 LogRegCoeff = LogRegCoeff,
                                                 DetectProb = 0.9,
                                                 verbose = T, showInfIndices = F,
                                                             LowerS = -1,
                                                             UpperS = 1),
                         
                       
                     grad = NULL,
                     ui = rbind(1, -1),
                     ci = c(a = -0.001, a = -0.001),
                     method = "Nelder-Mead")
cat("done!\n")

ML_ab <-  constrOptim(theta = c(a = ML_a$par[1], b = 0),
                     f = function(x) 
                       -AlleleFreqLogLikVar_4Par(FreqLengthDF = FreqLengthDF, 
                                                 L1TrueGivenEstList = L1TrueGivenEstList, 
                                                 L1MidPts = L1MidPts,
                                                 a = x[1], 
                                                 b = x[2], 
                                                 c = 0, 
                                                 d = 0, 
                                                 SD = NULL, 
                                                 N = PopSize, 
                                                 SampleSize = SampleSize,
                                                 blnIns = T, 
                                                 LogRegCoeff = LogRegCoeff,
                                                 DetectProb = 0.9,
                                                 verbose = T, showInfIndices = F,
                                                 LowerS = -1,
                                                 UpperS = 1),
                     
                     
                     grad = NULL,
                     ui = rbind(c(1, 0),  c(0, 1),   
                                c(-1, 0), c(0, -1)),
                     ci = c(a = -0.01, b = -10^-3, 
                            a = -0.01, b = -10^-3),
                     method = "Nelder-Mead")
cat("done!\n")

ML_abc <-  constrOptim(theta = c(a = ML_ab$par[1], b = ML_ab$par[2], 0),
                      f = function(x) 
                        -AlleleFreqLogLikVar_4Par(FreqLengthDF = FreqLengthDF, 
                                                  L1TrueGivenEstList = L1TrueGivenEstList, 
                                                  L1MidPts = L1MidPts,
                                                  a = x[1], 
                                                  b = x[2], 
                                                  c = x[3], 
                                                  d = 0, 
                                                  SD = NULL, 
                                                  N = PopSize, 
                                                  SampleSize = SampleSize,
                                                  blnIns = T, 
                                                  LogRegCoeff = LogRegCoeff,
                                                  DetectProb = 0.9,
                                                  verbose = T, showInfIndices = F,
                                                  LowerS = -1,
                                                  UpperS = 1),
                      
                      
                      grad = NULL,
                      ui = rbind(c(1, 0, 0),  c(0, 1, 0),  c(0, 0, 1), 
                                 c(-1, 0, 0), c(0, -1, 0), c(0, 0, -1)),
                      ci = c(a = -0.01, b = -10^-3, c = -0.01, 
                             a = -0.01, b = -10^-3, c = -0.01),
                      method = "Nelder-Mead")
cat("done!\n")


AlleleFreqLogLikVar_4Par(FreqLengthDF = FreqLengthDF, 
                         L1TrueGivenEstList = L1TrueGivenEstList, 
                         L1MidPts = L1MidPts,
                         a = -0.001, 
                         b = 0, 
                         c = 0, 
                         d = 0, 
                         SD = NULL, 
                         N = PopSize, 
                         SampleSize = SampleSize,
                         blnIns = T, 
                         LogRegCoeff = LogRegCoeff,
                         DetectProb = 0.9,
                         verbose = T, showInfIndices = F,
                         LowerS = -1,
                         UpperS = 1)

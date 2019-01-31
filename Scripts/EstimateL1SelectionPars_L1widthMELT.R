# The script below estimates the effect of L1 width on selection coefficients
# from 1000 genome data using insertion estimates obtained by MELT 
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
library(pracma)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

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
L1RefPath           <- 'D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv'
L1RefRangePath      <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
L1GRPath            <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
G1000SamplePath <- 'D:/L1polymORF/Data/1000GenomeSampleInfo.txt'
SelectTabOutPath    <- "D:/L1polymORF/Data/L1SelectionResults_L1widthMELT.csv"
SelectResultOutPath <- "D:/L1polymORF/Data/L1SelectionResults_L1widthMELT.RData"

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
MEInsCall$L1width <- sapply(MEInsCall$Info, GetLength)
#MEInsCall$SampleSize <- 1/min(MEInsCall$AF, na.rm = T)
MEInsCall$SampleSize <- 2 * MEInsSamplesize
MEInsCall$Freq <- MEInsCall$SampleSize * MEInsCall$AF
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
load(L1GRPath)
load(ChrLPath)
load(L1RefRangePath)
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

# Create a predictor variable for involvement in ectopic recombination
# L1Width      <- width(L1TotGR)
# hist(L1Width, breaks = seq(0, 6500, 100))
idxWidth <- which(!is.na(L1TotData$L1width))
L1TotData$RecPredict <- NA
L1TotData$RecPredict[idxWidth] <- 150000*sapply(L1TotData$L1width[idxWidth], function(x){
  idxMatch <- which.min(abs(x - DelVsL1Length$x))
  DelVsL1Length$y[idxMatch]
})

L1Width <- L1TotData$L1width[idxWidth]                             
L1Width[L1Width >= 4500] <- 4500
L1WidthOrder <- order(L1Width, decreasing = T)
OrderMatch   <- match(1:length(idxWidth), L1WidthOrder)
Deltas       <- c(L1Width[L1WidthOrder[-length(L1Width)]] - L1Width[L1WidthOrder[-1]],
                  L1Width[L1WidthOrder[length(L1Width)]])
DeltasSqProd <- 10^-7*Deltas^2 * (L1WidthOrder - 1)
Rev          <- length(idxWidth):1
L1WidthProd  <- cumsum(DeltasSqProd[Rev])[Rev]
L1TotData$RecPredict[idxWidth]  <- L1WidthProd[OrderMatch]
max(L1TotData$RecPredict[idxWidth])
plot(L1TotData$L1width, L1TotData$RecPredict)

# Number of L1 that are not fixed
Nnf <- nrow(L1TotData)

# Read information about 1000 genome samples
SampleInfo  <- read.table(G1000SamplePath, header = T)
SampleMatch <- match(SampleColumns, SampleInfo$sample)
Pops        <- SampleInfo$super_pop[SampleMatch]
NrS         <- length(SampleColumns)

# Define more genomic ranges
GeneGR <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
ExonGR <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
PromGR <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream = 10000)
CDSGR  <- cds(TxDb.Hsapiens.UCSC.hg19.knownGene)
IntronGRList   <- intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                      use.names = T)
FiveUTRGRList  <- fiveUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                       use.names = T)
ThreeUTRGRList <- threeUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                        use.names = T)
sum(width(GeneGR)/10^6) / sum(ChromLengthsHg19/10^6)

# Among overlapping genomic ranges, retain the longest
GeneGR <- UniqueGRanges(GeneGR)

cat("done!\n")


###################################################
#                                                 #
#   Fit effect of insertion length on selection   #
#                                                 #
###################################################

cat("\n********   Estimating effect of insertion length    **********\n")


# Create a matrix of predictor variables (L1 start and boolean variable for)
PredictMat <- L1TotData[, c("L1width", "blnFull", "RecPredict", 
                            "Freq", "SampleSize", "blnIns")]
                        
blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))

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
                          ci = c(a = -0.03, a = -0.03),
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
                          ci = c(a = -0.02, b = -10^(-6), 
                                 a = -0.02, b = -10^(-6)),
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
                      ci = c(a = -0.02, c = -10^(-3), 
                             a = -0.02, c = -10^(-3)),
                      method = "Nelder-Mead")
cat("done!\n")

# Determine maximum likelihood with 3 parameters (selection coefficient as 
# function of L1 start and indicator for full-length)
cat("Maximizing likelihood for three parameters ...")
ML_L1widthL1full <- constrOptim(theta = c(a = ML_L1width$par[1], b = ML_L1width$par[2], 
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
                  ci = c(a = -0.01, b = -10^(-6), c = -10^(-3), 
                         a = -0.02, b = -10^(-6), c = -10^(-3)),
                  method = "Nelder-Mead")
cat("done!\n")

# Determine maximum likelihood with 3 parameters (selection coefficient as 
# function of Recombination predictor and indicator for full-length)
# cat("Maximizing likelihood for three parameters ...")
# ML_L1RecL1full <- constrOptim(theta = c(a = ML_L1widthL1full$par[1],  
#                                           c = ML_L1widthL1full$par[3], d = ML_L1widthL1full$par[2]),
#                                 f = function(x) -AlleleFreqLogLik_4Par(
#                                   Freqs = round(L1TotData$Freq[!blnNA], 0),
#                                   Counts = rep(1, sum(!blnNA)),
#                                   Predict = PredictMat[!blnNA, 2:4],
#                                   a = x[1], b = 0, c = x[2], d = x[3], N = PopSize, 
#                                   SampleSize = L1TotData$SampleSize[!blnNA],
#                                   blnIns = L1TotData$blnIns[!blnNA], 
#                                   LogRegCoeff = LogRegL1Ref$coefficients,
#                                   DetectProb = L1TotData$DetectProb[!blnNA]),
#                                 grad = NULL,
#                                 ui = rbind(c(1, 0, 0),  c(0, 1, 0),  c(0, 0, 1), 
#                                            c(-1, 0, 0), c(0, -1, 0), c(0, 0, -1)),
#                                 ci = c(a = -0.01, b = -10^(-3), d = -10^(-6), 
#                                        a = -0.02, b = -10^(-3), d = -10^(-6)),
#                                 method = "Nelder-Mead")
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


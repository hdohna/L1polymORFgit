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
MeltInsPath         <- "D:/L1polymORF/Data/L1_SingleMELT_CombinedVcfs"
MeltDelPath         <- "D:/L1polymORF/Data/DEL.final_comp.vcf"
ChrLPath            <- 'D:/L1polymORF/Data/ChromLengthsHg19.Rdata'
InputPath           <- 'D:/L1polymORF/Data/SingletonAnalysis_unphased.RData'
L1RefPath           <- 'D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv'
L1RefRangePath      <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
RegrOutputPath      <- "D:/L1polymORF/Data/L1RegressionResults.RData"
SelectTabOutPath    <- "D:/L1polymORF/Data/L1SelectionResults_MELT_Single.csv"
SelectResultOutPath <- "D:/L1polymORF/Data/L1SelectionResults_MELT_Single.RData"

# Specify logistic rgeression coefficients for relationship between insert
# size and detection probability
L1SizeDetectCoeff <- c(a = 0.8, b = 10^-4)

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

# Coefficients for the probability to be reference, depending on L1 frequency
LogRegL1RefCoeff <- c(-4.706573, 9.737618)

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("\n\nLoading and processing data ...")

# Read in vcf file with MELT insertion calls
MEInsCallPerL1 <- read.table(MeltInsPath, as.is = T)

# Subset to get entries with genotype
idxWithGeno <- which(nchar(MEInsCallPerL1$Genotype) > 0)
MEInsCallPerL1 <- MEInsCallPerL1[idxWithGeno, ]

# Function to get numeric genotype and insertion length
GetGenoNum <- function(x){
  Split1 <- strsplit(x, "/")[[1]]
  Split2 <- strsplit(Split1[2], ":")[[1]][1]
  sum(as.numeric(c(Split1[1], Split2)))
}
GetLength <- function(x){
  Split1 <- strsplit(x, ";")[[1]]
  LengthPart <- grep("SVLEN=", Split1, value = T)
  if (length(LengthPart) > 0){
    as.numeric(strsplit(LengthPart, "=")[[1]][2])
  } else {
    NA
  }
}

# Add columns for numeric genotype and insertion length
MEInsCallPerL1$GenoNum <- sapply(MEInsCallPerL1$Genotype, GetGenoNum)
MEInsCallPerL1$L1width <- sapply(MEInsCallPerL1$INFO, GetLength)

# Get L1 ID and aggregate values per L1
MEInsCallPerL1$L1ID <- paste(MEInsCallPerL1$X.CHROM, MEInsCallPerL1$POS)
MEInsCall <- aggregate(MEInsCallPerL1$GenoNum, 
                       by = list(MEInsCallPerL1$X.CHROM, MEInsCallPerL1$POS), 
                       FUN = sum)
colnames(MEInsCall) <- c("CHROM", "POS", "Freq")
MEInsCall$L1ID <- paste(MEInsCall$CHROM, MEInsCall$POS)
L1IDmatch <- match(MEInsCall$L1ID, MEInsCallPerL1$L1ID)
MEInsCall$L1width <- MEInsCallPerL1$L1width[L1IDmatch]

# Add columns necessary for analysis 
MEInsCall$AF <- MEInsCall$Freq / MEInsSamplesize
MEInsCall$SampleSize <- 2 * MEInsSamplesize
MEInsCall$blnFull    <- MEInsCall$L1width >= MinLengthFullL1

# Create GRanges object for MEInsCall
MEInsCall$ChromName <- paste("chr", MEInsCall$CHROM, sep = "")
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
L1TotData$DetectProb[L1TotData$blnIns] <- 
  exp(L1SizeDetectCoeff[1] + 
        L1SizeDetectCoeff[2] * L1TotData$L1width[L1TotData$blnIns]) / 
  (1 + exp(L1SizeDetectCoeff[1] + 
             L1SizeDetectCoeff[2]*L1TotData$L1width[L1TotData$blnIns]))

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
                            LogRegCoeff = LogRegL1RefCoeff,
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
                            LogRegCoeff = LogRegL1RefCoeff,
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
                            LogRegCoeff = LogRegL1RefCoeff,
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
#                             ML_2Pars_L1count, 
                             ML_L1widthL1full, 
                             # ML_3Pars_L1countL1width,
                             # ML_3Pars_L1countL1full
#         ML_4Pars_L1countL1widthL1full
         ), function(x){
           c(AIC = GetAIC(x), Pars = GetParVals(x))
         }))
# Combine AIC values into one vector
AICTab <- cbind(data.frame(
            NrParameters = c(1, 
                             2, 
                             2, 
#                             2, 
                             3, 
                             # 3, 
                             # 3
#                             4
                             ),
            Predictor = c("none", 
                          "L1 width", 
                          "L1 full-length", 
#                          "L1count",
                          "L1 width and full-length",
                          # "L1 count and L1 start",
                          # "L1 count and L1 full"
  #                        "L1 start, L1 full-length, L1count"
                          ),
            stringsAsFactors = F),
            Cols2Append)
                     
# Save table with AIC
write.csv(AICTab, SelectTabOutPath)
save.image(SelectResultOutPath)

###################################################
#                                                 #
#   Fit effect of genic insertion on selection    #
#                                                 #
###################################################

# Create a matrix of predictor variables (L1 start and boolean variable for)
PredictMatGeneOL <- L1TotData[, c("blnOLExon", "blnOLIntron", "blnOLProm")]
PredictMatGeneOL2 <- L1TotData[, c("blnOLGene", "blnOLIntron", "blnOLProm")]
blnNA <- sapply(1:nrow(PredictMatGeneOL), function(x) any(is.na(PredictMatGeneOL[x,]))) |
  sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))

# Get maximum likelihood estimate for effect of exonic L1 on selection
cat("Estimate effect of exon overlap on selections ...")
ML_L1Exon <-  constrOptim(theta = c(a = ML_1Par$par, b = 0),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = round(L1TotData$Freq[!blnNA], 0),
                            Counts = rep(1, sum(!blnNA)),
                            Predict = PredictMatGeneOL[!blnNA,],
                            a = x[1], b = x[2], c = 0, d = 0, N = PopSize,
                            SampleSize = L1TotData$SampleSize[!blnNA],
                            blnIns = L1TotData$blnIns[!blnNA], 
                            LogRegCoeff = LogRegL1Ref$coefficients,
                            DetectProb = L1TotData$DetectProb[!blnNA]),
                          grad = NULL,
                           ui = rbind(c(1, 0),  c(0, 1),   
                                      c(-1, 0), c(0, -1)),
                           ci = c(a = -0.001, b = -10^(-2), 
                                  a = -0.001, b = -10^(-2)),
                           method = "Nelder-Mead")
cat("done!\n")

# Get maximum likelihood estimate for effect of intronic L1 on selection
cat("Estimate effect of intron overlap on selections ...")
ML_L1Intron <-  constrOptim(theta = c(a = ML_1Par$par, c = 0),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = round(L1TotData$Freq[!blnNA], 0),
                            Counts = rep(1, sum(!blnNA)),
                            Predict = PredictMatGeneOL[!blnNA,],
                            a = x[1], b = 0, c = x[2], d = 0, N = PopSize,
                            SampleSize = L1TotData$SampleSize[!blnNA],
                            blnIns = L1TotData$blnIns[!blnNA], 
                            LogRegCoeff = LogRegL1Ref$coefficients,
                            DetectProb = L1TotData$DetectProb[!blnNA]),
                          grad = NULL,
                          ui = rbind(c(1, 0),  c(0, 1),   
                                     c(-1, 0), c(0, -1)),
                          ci = c(a = -0.01, c = -10^(-2), 
                                 a = -0.01, c = -10^(-2)),
                          method = "Nelder-Mead")
cat("done!\n")

# Get maximum likelihood estimate for effect of intronic L1 on selection
cat("Estimate effect of promoter overlap on selections ...")
ML_L1Prom <-  constrOptim(theta = c(a = ML_1Par$par, d = 0),
                            f = function(x) -AlleleFreqLogLik_4Par(
                              Freqs = round(L1TotData$Freq[!blnNA], 0),
                              Counts = rep(1, sum(!blnNA)),
                              Predict = PredictMatGeneOL[!blnNA,],
                              a = x[1], b = 0, c = 0, d = x[2], N = PopSize,
                              SampleSize = L1TotData$SampleSize[!blnNA],
                              blnIns = L1TotData$blnIns[!blnNA], 
                              LogRegCoeff = LogRegL1Ref$coefficients,
                              DetectProb = L1TotData$DetectProb[!blnNA]),
                          grad = NULL,
                            ui = rbind(c(1, 0),  c(0, 1),   
                                       c(-1, 0), c(0, -1)),
                            ci = c(a = -0.01, c = -10^(-2), 
                                   a = -0.01, c = -10^(-2)),
                            method = "Nelder-Mead")
cat("done!\n")

# Get maximum likelihood estimate for effect of exonic L1 on selection
cat("Estimate effect of exon nad intron overlap on selections ...")
ML_L1ExonIntron <-  constrOptim(
                          theta = c(a = ML_1Par$par, 
                                    b = ML_L1Exon$par[2],
                                    c = ML_L1Intron$par[2]),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = round(L1TotData$Freq[!blnNA], 0),
                            Counts = rep(1, sum(!blnNA)),
                            Predict = PredictMatGeneOL[!blnNA,],
                            a = x[1], b = x[2], c = x[3], d = 0, N = PopSize,
                            SampleSize = L1TotData$SampleSize[!blnNA],
                            blnIns = L1TotData$blnIns[!blnNA], 
                            LogRegCoeff = LogRegL1Ref$coefficients,
                            DetectProb = L1TotData$DetectProb[!blnNA]),
                          grad = NULL,
                          ui = rbind(c(1, 0, 0),  c(0, 1, 0), c(0, 0, 1),   
                                     c(-1, 0, 0), c(0, -1, 0) , c(0, 0, -1)),
                          ci = c(a = -0.01, b = -10^(-2), c = -10^(-2), 
                                 a = -0.01, b = -10^(-2), c = -10^(-2)),
                          method = "Nelder-Mead")
cat("done!\n")

# Get maximum likelihood estimate for effect of exonic L1 on selection
cat("Estimate effect of exon nad intron overlap on selections ...")
ML_L1ExonIntronProm <-  constrOptim(
  theta = c(a = ML_L1ExonIntron$par[1], 
            b = ML_L1ExonIntron$par[2],
            c = ML_L1ExonIntron$par[3],
            d = 0),
f = function(x) -AlleleFreqLogLik_4Par(
  Freqs = round(L1TotData$Freq[!blnNA], 0),
  Counts = rep(1, sum(!blnNA)),
  Predict = PredictMatGeneOL[!blnNA,],
  a = x[1], b = x[2], c = x[3], d = x[4], N = PopSize,
  SampleSize = L1TotData$SampleSize[!blnNA],
  blnIns = L1TotData$blnIns[!blnNA], 
  LogRegCoeff = LogRegL1Ref$coefficients,
  DetectProb = L1TotData$DetectProb[!blnNA]),
grad = NULL,
  ui = rbind(c(1, 0, 0, 0),  c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1),  
             c(-1, 0, 0, 0), c(0, -1, 0, 0) , c(0, 0, -1, 0), c(0, 0, 0, -1)),
  ci = c(a = -0.01, b = -10^(-2), c = -10^(-2), d = -10^(-2), 
         a = -0.01, b = -10^(-2), c = -10^(-2), d = -10^(-2)),
  method = "Nelder-Mead")
cat("done!\n")

# Get maximum likelihood estimate for effect of exonic L1 on selection
cat("Estimate effect of exon or intron overlap on selections ...")
ML_L1PromOrIntron <-  constrOptim(theta = c(a = ML_1Par$par, 
                                            b = ML_L1Exon$par[2],
                                            c = ML_L1Intron$par[2]),
                                  f = function(x) -AlleleFreqLogLik_4Par(
                                    Freqs = round(L1TotData$Freq[!blnNA], 0),
                                    Counts = rep(1, sum(!blnNA)),
                                    Predict = PredictMatGeneOL[!blnNA,],
                                    a = x[1], b = x[2], c = x[3], d = x[3], N = PopSize,
                                    SampleSize = L1TotData$SampleSize[!blnNA],
                                    blnIns = L1TotData$blnIns[!blnNA], 
                                    LogRegCoeff = LogRegL1Ref$coefficients,
                                    DetectProb = L1TotData$DetectProb[!blnNA]),
                                  grad = NULL,
                                  ui = rbind(c(1, 0, 0),  c(0, 1, 0), c(0, 0, 1),   
                                             c(-1, 0, 0), c(0, -1, 0) , c(0, 0, -1)),
                                  ci = c(a = -0.01, b = -10^(-2), c = -10^(-2), 
                                         a = -0.01, b = -10^(-2), c = -10^(-2)),
                                  method = "Nelder-Mead")
cat("done!\n")


# Get columns of AIC and parameter values
Cols2Append <- t(sapply(list(ML_1Par, ML_L1Exon, ML_L1Intron, ML_L1Prom, 
                             ML_L1ExonIntron,
                             ML_L1ExonIntronProm,
                             ML_L1PromOrIntron), 
                        function(x){
                               c(NrParameters = GetNPar(x), AIC = GetAIC(x), 
                                 Pars = GetParVals(x))
                             }))
# Combine AIC values into one vector
AICTabGene <- cbind(data.frame(
  Predictor = c("none", "Exon", "Intron", "Promoter",
                "Exon and intron", 
                "Exon, intron, and promoter",
                "Exon, intron or promoter"),
  stringsAsFactors = F),
  Cols2Append)

# Save table with AIC
write.csv(AICTabGene, SelectGenTabOutPath)

###################################################
#                                                 #
#   Plot density vs. selection coefficient        #
#                                                 #
###################################################

# Create a vector of selection coefficients
SCoeffVect <- c(Promoter = ML_L1ExonIntron$par[1],
                Exon = sum(ML_L1ExonIntron$par[c(1, 2)]),
                Intron = sum(ML_L1ExonIntron$par[c(1, 3)]),
                Intergenic = ML_L1ExonIntron$par[1])
names(SCoeffVect) <- sapply(names(SCoeffVect), 
                            function(x) strsplit(x, "\\.")[[1]][1])

# Plot selection coefficient against 
if (!all(names(SCoeffVect) == colnames(InsPerbp))){
  stop("Selection coefficients and L1 densities are not in same order!")
}
if (!all(names(SCoeffVect) == names(MeanFreqs))){
  stop("Selection coefficients and L1 frequencies are not in same order!")
}

# Get sample size and create a range of s-values
SSize <- 2 * MEInsSamplesize
SVals  <- seq(-0.0025, -0.00001, 0.00001)

# Plot probability for inclusion versus number of LINE-1 per Mb
ProbL1 <- sapply(SVals, function(x) ProbAlleleIncluded(x,N = PopSize, SampleSize = 2*2504))
par(oma = c(7, 1, 0, 2), mfrow = c(2, 1), mai = c(0.5, 1, 0.5, 1))
plot(SCoeffVect, InsPerbp[2,], ylab = "LINE-1s per Mb", 
     xlab = "", ylim = c(0, 3), xlim = c(-0.0025, 0), main = "A")
text(SCoeffVect, InsPerbp[2,] + 2*10^(-1), names(SCoeffVect))
par(new = TRUE)
plot(SVals, ProbL1, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side = 4)
mtext("Inclusion probability", 4, line = 3)

# Plot expected frequency versus observed mean frequency
ExpL1 <- sapply(SVals, function(x) ExpAlleleFreq(x, N = PopSize, SampleSize = 2*2504))
plot(SCoeffVect, MeanFreqs*SSize, ylab = "Mean LINE-1 frequency", 
     xlab = "", xlim = c(-0.0025, 0.0001), main = "B")
text(SCoeffVect + c(0.0002, 0, -0.0001, -0.0002), MeanFreqs*SSize + 10, names(SCoeffVect))
lines(SVals, ExpL1)
mtext("Selection coefficient", 1, line = 3)
CreateDisplayPdf('D:/L1polymORF/Figures/SelectionPerRegion_MELT.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7, width = 7)

###################################################
#                                                 #
#   Plot frequency vs. insertion length           #
#                                                 #
###################################################

# # Create a vector of L1 start classes
# L1TotData$L1widthClass <- cut(L1TotData$L1width, breaks = 
#                                   seq(0, 7000, 1000))
# MEInsCall$L1widthClass <- cut(MEInsCall$L1width, breaks = 
#                                    seq(0, 7000, 1000))
# 
# MEInsCall$Freq
# # Get mean L1 frequency per start
# L1widthAggregated <- aggregate(L1TotData[,c("L1width", "L1Freq")], 
#                                by = list(L1TotData$L1widthClass), 
#                                FUN = function(x) mean(x, na.rm = T))
# L1widthAggregated_Ins <- aggregate(MEInsCall[,c("L1width", "AF")], 
#                                by = list(MEInsCall$L1widthClass), 
#                                FUN = function(x) mean(x, na.rm = T))
# plot(L1widthAggregated_Ins$L1width, L1widthAggregated_Ins$AF)
# 
# # Get sample size and create a range of s-values
# SSize <- 2 * MEInsSamplesize
# StartVals  <- seq(0, 6000, 100)
# Full       <- StartVals == 6000
# SVals <- ML_L1widthL1full$par[1] + ML_L1widthL1full$par[2]*StartVals +
#   ML_L1widthL1full$par[3]*Full
# 
# # Plot expected frequency versus observed mean frequency
# ExpL1width <- sapply(SVals, function(x) ExpAlleleFreq(x, N = PopSize, 
#                                                       SampleSize = 2*MEInsSamplesize))
# par( mfrow = c(1, 1))
# plot(L1widthAggregated$L1width, 
#      L1widthAggregated$L1Freq, xlab = "LINE-1 length",
#      ylab = "Mean LINE-1 frequency")
# lines(StartVals, ExpL1width )
# mtext("Selection coefficient", 1, line = 3)
# CreateDisplayPdf('D:/L1polymORF/Figures/FreqVsL1width_MELT.pdf',
#                  PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
#                  height = 7, width = 7)

###################################################
#                                                 #
#   Fit effect of strandedness on selection       #
#                                                 #
###################################################

# Create a matrix of predictor variables (L1 start and boolean variable for)
PredictMatWithinGene <- L1TotData[L1TotData$blnOLGene & !blnNA , 
                                 c( "blnOLGeneSameStrand", "blnOLGene", "blnOLGene")]

# Estimate maximum likelihood for a single selection coefficient
sum(L1TotData$blnOLGene)
colSums(PredictMatWithinGene)
ML_1Par_gene <- constrOptim(theta = c(a = 0),
                            f = function(x) -AlleleFreqLogLik_4Par(
                              Freqs = round(L1TotData$Freq[L1TotData$blnOLGene & !blnNA], 0),
                              Counts = rep(1, sum(L1TotData$blnOLGene & !blnNA)),
                              Predict = PredictMatWithinGene,
                              a = x[1], b = 0, c = 0, d = 0, N = PopSize,
                              SampleSize = L1TotData$SampleSize[L1TotData$blnOLGene  & !blnNA ],
                              blnIns = L1TotData$blnIns[L1TotData$blnOLGene & !blnNA], 
                              LogRegCoeff = LogRegL1Ref$coefficients,
                              DetectProb = L1TotData$DetectProb[L1TotData$blnOLGene & !blnNA]),
                            grad = NULL,
                            ui = rbind(1,-1),
                            ci = c(a = -0.001, a = -0.001),
                            method = "Nelder-Mead")

# Get maximum likelihood estimate for effect of exonic L1 on selection
cat("Estimate effect of same strand overlap on selections ...")
ML_L1SameStrand <-  constrOptim(theta = c(a = ML_1Par_gene$par, b = 0),
                                f = function(x) -AlleleFreqLogLik_4Par(
                                  Freqs = round(L1TotData$Freq[L1TotData$blnOLGene & !blnNA], 0),
                                Counts = rep(1, sum(L1TotData$blnOLGene & !blnNA)),
                                Predict = PredictMatWithinGene,
                                a = x[1], b = x[2], c = 0, d = 0, N = PopSize,
                                SampleSize = L1TotData$SampleSize[L1TotData$blnOLGene  & !blnNA ],
                                blnIns = L1TotData$blnIns[L1TotData$blnOLGene & !blnNA], 
                                LogRegCoeff = LogRegL1Ref$coefficients,
                                DetectProb = L1TotData$DetectProb[L1TotData$blnOLGene & !blnNA]),
                          grad = NULL,
                          ui = rbind(c(1, 0),  c(0, 1),   
                                     c(-1, 0), c(0, -1)),
                          ci = c(a = -0.01, b = -10^(-2), 
                                 a = -0.01, b = -10^(-2)),
                          method = "Nelder-Mead")
cat("done!\n")

# Get columns of AIC and parameter values
Cols2Append <- t(sapply(list(ML_1Par_gene, ML_L1SameStrand), 
                        function(x){
                          c(NrParameters = GetNPar(x), AIC = GetAIC(x), 
                            Pars = GetParVals(x))
                        }))
# Combine AIC values into one vector
AICTabWithinGene <- cbind(data.frame(
  Predictor = c("none", "SameStrand"),
  stringsAsFactors = F),
  Cols2Append)

# Save table with AIC
write.csv(AICTabWithinGene, SelectWithinGenTabOutPath)

###################################################
#                                                 #
#   Fit effect of singleton coef. on selection    #
#                                                 #
###################################################

# Create a matrix of predictor variables 
PredictMat <- L1SingletonCoeffs[, c("coef", "coef", "coef")]
blnNA <- sapply(1:nrow(L1SingletonCoeffs), function(x) any(is.na(PredictMat[x,])))

# Determine maximum likelihood with one parameter (selection coefficient)
cat("Maximizing likelihood for one parameter (selection coefficient) ...")
ML_1Par_coef <- constrOptim(
  theta = c(a = ML_1Par$par),
  f = function(x) -AlleleFreqLogLik_4Par(
    Freqs = (L1SingletonCoeffs$Freq * 2*2504)[!blnNA], 
    Counts = rep(1, sum(!blnNA)), 
    Predict = PredictMat[!blnNA,], 
    a = x[1], b = 0, c = 0, d = 0, N = PopSize,
    SampleSize = rep(2*2504, sum(!blnNA)),
    blnIns = rep(T, sum(!blnNA)), 
    LogRegCoeff = LogRegL1Ref$coefficients,
    DetectProb = rep(0.9, sum(!blnNA))),
  grad = NULL,
  ui = rbind(1,-1),
  ci = c(a = -0.03, a = -0.03),
  method = "Nelder-Mead")

# Determine maximum likelihood with an intercept and one parameter for thr 
# selection coefficient
ML_2Pars_L1coef <- constrOptim(
  theta = c(a = ML_1Par_coef$par, b = 0),
  f = function(x) -AlleleFreqLogLik_4Par(
    Freqs = (L1SingletonCoeffs$Freq * 2*2504)[!blnNA],
    Counts = rep(1, sum(!blnNA)),
    Predict = PredictMat[!blnNA,], 
    a = x[1], b = x[2], c = 0, d = 0, N = PopSize,
    SampleSize = rep(2*2504, sum(!blnNA)),
    blnIns = rep(T, sum(!blnNA)), 
    LogRegCoeff = LogRegL1Ref$coefficients,
    DetectProb = rep(0.9, sum(!blnNA))),
  grad = NULL,
  ui = rbind(c(1, 0),  c(0, 1),  
             c(-1, 0), c(0, -1)),
  ci = c(a = -0.01, b = -2*10^(-3), 
         a = -0.01, b = -2*10^(-3)),
  method = "Nelder-Mead")
cat("done!\n")

# Get columns of AIC and parameter values
Cols2Append <- t(sapply(list(ML_1Par_coef, ML_2Pars_L1coef), 
                        function(x){
                          c(NrParameters = GetNPar(x), AIC = GetAIC(x), 
                            Pars = GetParVals(x))
                        }))
# Combine AIC values into one vector
AICTabSingleton <- cbind(data.frame(
  Predictor = c("none", "Signleton coefficient"),
  stringsAsFactors = F),
  Cols2Append)

# Save table with AIC
write.csv(AICTabSingleton, SelectSingletonTabOutPath)

# Save everything
save.image(SelectResultOutPath)

# The script below estimates selection coefficients of L1 from the 
# 1000 genome data using approximate bayesian computing (ABC) 

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
library(abc)

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

# Vector of frequency counts
FreqCountV <- 1:30

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("\n\nLoading and processing data ...")

# Load simulation results
load("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1Simulated_AdditionalInfo_MELT.RData")

# Source start script again
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

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

# Group L1 width into bins
MEInsCall$L1LengthBins <- cut(MEInsCall$L1width, breaks = seq(0, 6500, 500))

# Get mean L1 frequency per length
L1WidthAggregated <- AggDataFrame(MEInsCall, 
                                  GroupCol = "L1LengthBins", 
                                  MeanCols = c("L1width", "Freq", "blnFull"), 
                                  LengthCols = "Freq",
                                  VarCols = "Freq")

# Get per L1 width bin the site-frequency spectrum (SFS) and calculate a 
# multinomial probability per SFS
LBinsUnique <- unique(MEInsCall$L1LengthBins)
LBinsUnique <- LBinsUnique[!is.na(LBinsUnique)]
SFS_list <- lapply(LBinsUnique, function(x){
  blnX   <- MEInsCall$L1LengthBins == x
  Counts <- MEInsCall$Freq[blnX]
  Counts <- Counts[!is.na(Counts)]
  if (length(Counts) > 0){
    SFS <- sapply(FreqCountV, function(y) sum(Counts == y))
    SFS <- c(SFS, sum(Counts > max(Counts)))
  }
})
names(SFS_list) <- LBinsUnique

# Match data fram to SFS_list
binMatch <- match(names(SFS_list), L1WidthAggregated$L1LengthBins)
L1WidthAggregated <- L1WidthAggregated[binMatch,]

# Create a vector of observed multinomial frequency probabilities
DMulti <- sapply(SFS_list, function(x) log(dmultinom(x, sum(x), x)))

# Create matrix of a combination of parameter values
a <- seq(-10^-2, 10^-2, 5*10^-3)
b <- seq(-10^-5, 10^-5, 5*10^-6)
c <- seq(-10^-2, 10^-2, 5*10^-3)
ParCombos <- data.frame(aVals = rep(rep(a, length(b)), length(c)),
           bVals = rep(rep(b, each = length(a)), length(c)),
           cVals = rep(c, each = length(a) * length(b)))

# Calculate vector of multinomial probabilities for each parameter
# combination
x <- 1
sPredictors <- as.matrix(cbind(1, L1WidthAggregated[,c("L1width_mean", "blnFull_mean")]))
SimMat <- t(sapply(1:20, function(x){
  Pars <- t(as.vector(ParCombos[x,]))
  sVals <- sPredictors %*% Pars
  sapply(1:nrow(sPredictors), function(i){
    SimulateMultinomSFSSelection(SFS_list[[i]], FreqCountV, s = sVals[i], 
                                 SampleSize = 5000,
                                  N = 10^4, NAdd = 100, 
                                  NGenPerBatch = 10^2, MaxNGen = 10^4,
                                  PDiffThresh = 0.3)
  })
}))
blnInf <- apply(SimMat, 1, FUN = function(x) any(is.infinite(x)))
dim(SimMat)
abc(DMulti, ParCombos[which(!blnInf), ], SimMat[!blnInf, ], tol = 0.01, method = "loclinear")

###################################################
#                                                 #
#   Fit effect of insertion length on selection   #
#                                                 #
###################################################


# Save everything
save.image(SelectResultOutPath)

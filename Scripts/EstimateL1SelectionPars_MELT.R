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
load(InputPath)
load(L1GRPath)
load(ChrLPath)
load(L1RefRangePath)
load(RegrOutputPath)

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

# Make genomic ranges for L1SingletonCoeffs
L1SingletonCoeffs$chromosome <- paste("chr", L1SingletonCoeffs$Chrom, sep = "")
L1SingletonCoeffs_GR <- makeGRangesFromDataFrame(L1SingletonCoeffs, 
                                                 seqnames.field = "chromosome",
                                                 start.field = "Pos",
                                                 end.field = "Pos")

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

##########################################
#                                        #
#  Add columns to L1SingletonCoeffs      #
#                                        #
##########################################

cat("Add columns to L1SingletonCoeffs ...")

# Turn factors into numeric values
L1SingletonCoeffs$L1Start <- as.numeric(as.character(L1SingletonCoeffs$L1Start))
L1SingletonCoeffs$L1End <- as.numeric(as.character(L1SingletonCoeffs$L1End))

# Indicator for full-length
L1SingletonCoeffs$blnFull <- L1SingletonCoeffs$L1Start <= 3 &
  L1SingletonCoeffs$L1End >= MinLengthFullL1
sum(L1SingletonCoeffs$InsLength <= 100)

# Indicator for significant effect
L1SingletonCoeffs$blnSig <- p.adjust(L1SingletonCoeffs$Pr...z..) < FDR
hist(L1SingletonCoeffs$Pr...z.., breaks = seq(0, 1, 0.005))

# Indicator for positive selection
L1SingletonCoeffs$blnSelect <- L1SingletonCoeffs$blnSig &
  L1SingletonCoeffs$coef < 0

# Indicator for negative selection
L1SingletonCoeffs$blnNegSelect <- L1SingletonCoeffs$blnSig &
  L1SingletonCoeffs$coef > 0
sum(L1SingletonCoeffs$blnNegSelect)
sum(L1SingletonCoeffs$blnSelect)

# Indicator fo selection (+1 = positive, -1 = negative, 0 = neutral)
L1SingletonCoeffs$SelectInd <- 0
L1SingletonCoeffs$SelectInd[L1SingletonCoeffs$blnSelect]     <- 1
L1SingletonCoeffs$SelectInd[L1SingletonCoeffs$blnNegSelect]  <- -1

# Caclulate distance to genes
L1SingletonCoeffs$Dist2Gene <- Dist2Closest(L1SingletonCoeffs_GR, GeneGR)
L1SingletonCoeffs$blnOLGene <- L1SingletonCoeffs$Dist2Gene == 0

# Caclulate logarithm of distance to genes
L1SingletonCoeffs$LogDist2Gene <- log(L1SingletonCoeffs$Dist2Gene + 0.1)

# Standardize SE ratio to one
L1SingletonCoeffs$SE_RatioSt <- 1/L1SingletonCoeffs$se.coef./
  mean(1/L1SingletonCoeffs$se.coef.)

# Add boolean indicators for overlap
L1SingletonCoeffs$blnOLGene  <- overlapsAny(L1SingletonCoeffs_GR, GeneGR, ignore.strand = T)
L1SingletonCoeffs$blnOLProm     <- overlapsAny(L1SingletonCoeffs_GR, PromGR, ignore.strand = T)
L1SingletonCoeffs$blnOLExon     <- overlapsAny(L1SingletonCoeffs_GR, ExonGR, ignore.strand = T)
L1SingletonCoeffs$blnOLIntron   <- L1SingletonCoeffs$blnOLGene & (!L1SingletonCoeffs$blnOLExon)
L1SingletonCoeffs$blnOLIntergen <- !(L1SingletonCoeffs$blnOLGene | L1SingletonCoeffs$blnOLProm)
L1SingletonCoeffs$L1StartNum <- as.numeric(as.character(L1SingletonCoeffs$L1Start))
L1SingletonCoeffs$L1EndNum   <- as.numeric(as.character(L1SingletonCoeffs$L1End))
L1SingletonCoeffs$blnFull    <- L1SingletonCoeffs$L1StartNum <= 1 & L1SingletonCoeffs$L1EndNum >= 6000

# Add info about overlapping genes
L1coeff_Gene_OL <- findOverlaps(L1SingletonCoeffs_GR, GeneGR, ignore.strand = T)
L1SingletonCoeffs$idxGene   <- NA
L1SingletonCoeffs$GeneWidth <- NA
L1SingletonCoeffs$GeneID <- NA
L1SingletonCoeffs$blnOLGeneSameStrand <- NA
L1SingletonCoeffs$idxGene[L1coeff_Gene_OL@from] <- L1coeff_Gene_OL@to
L1SingletonCoeffs$GeneWidth[L1coeff_Gene_OL@from] <- width(GeneGR)[L1coeff_Gene_OL@to]
L1SingletonCoeffs$GeneID[L1coeff_Gene_OL@from] <- 
  GeneGR@elementMetadata@listData$gene_id[L1coeff_Gene_OL@to]
L1SingletonCoeffs$blnOLGeneSameStrand[L1coeff_Gene_OL@from] <- 
  L1SingletonCoeffs$L1Strand[L1coeff_Gene_OL@from] == as.vector(strand(GeneGR))[L1coeff_Gene_OL@to]

# Check out properties of L1 with signal of positive selectiom
L1SingletonCoeffs[L1SingletonCoeffs$blnSelect,]
fisher.test(L1SingletonCoeffs$blnSelect, (L1SingletonCoeffs$blnOLIntron & 
                                            L1SingletonCoeffs$blnOLGeneSameStrand))
mean((L1SingletonCoeffs$blnOLIntron & 
        L1SingletonCoeffs$blnOLGeneSameStrand))

# Standardize selection coefficients
CoeffAggMean <- aggregate(coef ~ Freq, data = L1SingletonCoeffs, FUN = mean)
CoeffAggVar  <- aggregate(coef ~ Freq, data = L1SingletonCoeffs, FUN = var)
CoeffAggN    <- aggregate(coef ~ Freq, data = L1SingletonCoeffs, FUN = length)
CoeffAggMerge <- merge(CoeffAggMean, CoeffAggVar, by = 'Freq')
CoeffAggMerge <- merge(CoeffAggMerge, CoeffAggN, by = 'Freq')
colnames(CoeffAggMerge)[2:4] <- c("Mean", "Var", "N")
CoeffAggMerge$StDev <- sqrt(CoeffAggMerge$Var)

FreqMatch <- match(L1SingletonCoeffs$Freq, CoeffAggMerge$Freq)
L1SingletonCoeffs$MeanCof <- CoeffAggMerge$Mean[FreqMatch]
L1SingletonCoeffs$StDevCof <- CoeffAggMerge$StDev[FreqMatch]
L1SingletonCoeffs$CoefSt <- (L1SingletonCoeffs$coef - L1SingletonCoeffs$MeanCof) / 
  L1SingletonCoeffs$StDevCof
cat("done!\n")

####################################################
#                                                  #
#   Overview of L1 intersection with features      #
#                                                  #
####################################################

# Indicator variable for intersection with various GRanges
L1TotData$blnOLGene  <- overlapsAny(L1TotGR, GeneGR, ignore.strand = T)
L1TotData$blnOLGeneSameStrand <- overlapsAny(L1TotGR, GeneGR)
L1TotData$blnOLProm     <- overlapsAny(L1TotGR, PromGR, ignore.strand = T)
L1TotData$blnOLExon     <- overlapsAny(L1TotGR, ExonGR, ignore.strand = T)
L1TotData$blnOLIntron   <- L1TotData$blnOLGene & (!L1TotData$blnOLExon)
L1TotData$blnOLIntergen <- !(L1TotData$blnOLGene | L1TotData$blnOLProm)
  
# Create a variable indicating insertion type
L1TotData$InsType                          <- "Intergenic"
L1TotData$InsType[L1TotData$blnOLProm]     <- "Promoter"
  L1TotData$InsType[L1TotData$blnOLExon]   <- "Exon"
  L1TotData$InsType[L1TotData$blnOLIntron] <- "Intron"

# Perform pairwise Wilcoxon test for differences in L1 frequencies
pairwise.wilcox.test(L1TotData$L1Freq, L1TotData$InsType,
                     p.adjust.method = "BH")
# Average mean frequency
MeanFreqAgg <- aggregate(L1Freq ~ InsType, data = L1TotData, FUN = mean)
VarFreqAgg  <- aggregate(L1Freq ~ InsType, data = L1TotData, FUN = var)
L1TotData$Dummy <- 1
NAgg  <- aggregate(Dummy ~ InsType, data = L1TotData, FUN = sum)
StErr <- sqrt(VarFreqAgg$Frequency / NAgg$Dummy)

# Indicator variable for intersection with reference L1
blnOLGene_RefL1  <- overlapsAny(L1GRanges, GeneGR, ignore.strand = T)
blnOLGeneSameStrand_RefL1 <- overlapsAny(L1GRanges, GeneGR)
blnOLProm_RefL1   <- overlapsAny(L1GRanges, PromGR, ignore.strand = T)
blnOLExon_RefL1   <- overlapsAny(L1GRanges, ExonGR, ignore.strand = T)
blnOLIntron_RefL1 <- blnOLGene_RefL1 & (!blnOLExon_RefL1)

# Get number of insertions per bp
GeneTot     <- sum(width(GeneGR))
ExonTot     <- sum(width(ExonGR))
IntronTot   <- GeneTot - ExonTot
PromTot     <- sum(width(PromGR))
IntergenTot <- sum(as.numeric(ChromLengthsHg19)) - GeneTot - PromTot #- EnhancerTot

# Get mean frequency of L1 in different functional regions
MeanFreqs <- c(
               Promoter = mean(L1TotData$L1Freq[L1TotData$blnOLProm], na.rm = T),
               Exon = mean(L1TotData$L1Freq[L1TotData$blnOLExon], na.rm = T),
               Intron = mean(L1TotData$L1Freq[L1TotData$blnOLIntron], na.rm = T),
               Intergenic = mean(L1TotData$L1Freq[L1TotData$blnOLIntergen], na.rm = T)
               )

# Plot distn of frequency of L1 in different functional regions
par(mfrow = c(1, 1))
hist(L1TotData$L1Freq[L1TotData$blnOLProm], breaks = seq(0, 1, 0.01))
hist(L1TotData$L1Freq[L1TotData$blnOLExon], breaks = seq(0, 1, 0.01))
hist(L1TotData$L1Freq[L1TotData$blnOLIntron], breaks = seq(0, 1, 0.01))
hist(L1TotData$L1Freq[L1TotData$blnOLIntergen], breaks = seq(0, 1, 0.01))
hist(sqrt(-log10(L1TotData$L1Freq[L1TotData$blnOLProm])))
hist(-log10(L1TotData$L1Freq[L1TotData$blnOLExon]))
hist(log10(L1TotData$L1Freq[L1TotData$blnOLIntron]))
hist(log10(L1TotData$L1Freq[L1TotData$blnOLIntergen]))

# Get number of L1 per Mb in different functional regions
InsPerbp <- 10^6 * rbind(
  c(
    Promoter = sum(blnOLProm_RefL1) / PromTot,
    Exon = sum(blnOLExon_RefL1) / ExonTot,
    Intron = sum(blnOLIntron_RefL1) / IntronTot,
    Intergenic = sum(!(blnOLGene_RefL1 | blnOLProm_RefL1)) / IntergenTot
    ),
  c(
    Promoter = sum(L1TotData$blnOLProm) / PromTot,
    Exon = sum(L1TotData$blnOLExon) / ExonTot,
    Intron = sum(L1TotData$blnOLIntron) / IntronTot,
    Intergenic = sum(!(L1TotData$blnOLGene | L1TotData$blnOLProm)) / IntergenTot
    )
  
)
InsPerbp[1,] / InsPerbp[2,]

###################################################
#                                                 #
#   Fit effect of insertion length on selection   #
#                                                 #
###################################################

cat("\n********   Estimating effect of insertion length    **********\n")

# Match summary ranges to L1 ranges of 1000 genome data
L1SummaryOL <- findOverlaps(L1TotGR, SummaryGR)
all(L1SummaryOL@from %in% 1:nrow(L1TotData))
blnNoDupl <- !duplicated(L1SummaryOL@from)
L1TotData$L1Count <- NA
L1TotData$L1Count[L1SummaryOL@from[blnNoDupl]] <- 
  DataPerSummaryGR$L1Count[L1SummaryOL@to[blnNoDupl]]

# Get distance to nearest other L1
Dist2Nearest <- distanceToNearest(L1TotGR)
L1TotData$Dist2Nearest <- Dist2Nearest@elementMetadata@listData$distance
max(L1TotData$Dist2Nearest)

# Create a matrix of predictor variables (L1 start and boolean variable for)
PredictMat <- L1TotData[, c("L1Count", "L1width", "blnFull", "Freq", "SampleSize",
                            "blnIns")]
blnNA <- sapply(1:nrow(L1TotData), function(x) any(is.na(PredictMat[x,])))
sum(!blnNA)
max(L1TotData$Freq / L1TotData$SampleSize, na.rm = T)

# Plot log-likelihood for different selection coefficients
# aVals <- seq(-0.0021, 0.003, 0.0001)
# LikVals <- sapply(aVals, function(x) {
#   print(x)
#   LL_FPrime = AlleleFreqLogLik_4Par(
#     Freqs = round(L1TotData$Freq[!blnNA], 0),
#     Counts = rep(1, sum(!blnNA)),
#     Predict = PredictMat[!blnNA, 1:3],
#     a = x, b = 0, c = 0, d = 0, N = 10^4,
#     SampleSize = L1TotData$SampleSize[!blnNA], 
#     blnIns = L1TotData$blnIns[!blnNA], 
#     DetectProb = 0.9)
#   })
# par(mfrow = c(1, 1))
# plot(aVals, LikVals, type = "l", col = "red")
# plot(aVals, LikVals, type = "l", col = "red", 
#      xlim = c(-0.0005, 0), ylim)

# Estimate maximum likelihood for a single selection coefficient
cat("Estimate maximum likelihood for a single selection coefficient\n")
ML_1Par <-  constrOptim(theta = c(a = 0),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = round(L1TotData$Freq[!blnNA], 0),
                            Counts = rep(1, sum(!blnNA)),
                            Predict = PredictMat[!blnNA, 1:3],
                            a = x[1], b = 0, c = 0, d = 0, N = 10^4,
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
ML_L1width <-  constrOptim(theta = c(a = ML_1Par$par, c = 0),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = round(L1TotData$Freq[!blnNA], 0),
                            Counts = rep(1, sum(!blnNA)),
                            Predict = PredictMat[!blnNA, 1:3],
                            a = x[1], b = 0, c = x[2], d = 0, N = 10^4, 
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
ML_L1full <-  constrOptim(theta = c(a = ML_1Par$par, d = 0),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = round(L1TotData$Freq[!blnNA], 0),
                            Counts = rep(1, sum(!blnNA)),
                            Predict = PredictMat[!blnNA, 1:3],
                            a = x[1], b = 0, c = 0, d = x[2], N = 10^4, 
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
ML_L1widthL1full <- constrOptim(theta = c(a = ML_L1width$par[1], b = ML_L1width$par[2], 
                                          c = ML_L1full$par[2]),
                  f = function(x) -AlleleFreqLogLik_4Par(
                    Freqs = round(L1TotData$Freq[!blnNA], 0),
                    Counts = rep(1, sum(!blnNA)),
                    Predict = PredictMat[!blnNA, 1:3],
                    a = x[1], b = 0, c = x[2], d = x[3], N = 10^4, 
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
#  Fit effect of L1 density on selection          #
#                                                 #
###################################################


# Determine maximum likelihood with 3 parameters (selection coefficient as 
# function of L1 start and indicator for full-length)
cat("Maximizing likelihood for L1 count ...")
ML_2Pars_L1count <- constrOptim(
     theta = c(a = ML_1Par$par, b = 0),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = round(L1TotData$Freq[!blnNA], 0),
                            Counts = rep(1, sum(!blnNA)),
                            Predict = PredictMat[!blnNA, 1:3],
                            a = x[1], b = x[2], c = 0, d = 0, N = 10^4, 
                            SampleSize = L1TotData$SampleSize[!blnNA],
                            blnIns = L1TotData$blnIns[!blnNA], 
                            LogRegCoeff = LogRegL1Ref$coefficients,
                            DetectProb = L1TotData$DetectProb[!blnNA]),
                          grad = NULL,
                          ui = rbind(c(1, 0),  c(0, 1),  
                                     c(-1, 0), c(0, -1)),
     # ci = c(a = -0.01, b = -10^(-3), 
     #        a = -0.01, b = -2*10^(-3)),
     ci = c(a = -0.01, b = -10^(-9), 
            a = -0.01, b = -10^(-9)),
     method = "Nelder-Mead")

# Maximum likelihood estimate for effect of L1 density and full-length L1
ML_3Pars_L1countL1full <- constrOptim(
  theta = c(a = ML_2Pars_L1count$par[1], ML_2Pars_L1count$par[2], 
            d = ML_L1full$par[2]), 
  f = function(x) -AlleleFreqLogLik_4Par(
    Freqs = round(L1TotData$Freq[!blnNA], 0),
    Counts = rep(1, sum(!blnNA)),
    Predict = PredictMat[!blnNA, 1:3],
    a = x[1], b = x[2], c = 0, d = x[3], N = 10^4, 
    SampleSize = L1TotData$SampleSize[!blnNA],
    blnIns = L1TotData$blnIns[!blnNA], 
    LogRegCoeff = LogRegL1Ref$coefficients,
    DetectProb = L1TotData$DetectProb[!blnNA]),
  grad = NULL,
  ui = rbind(c(1, 0, 0),  c(0, 1, 0), c(0, 0, 1),     
             c(-1, 0, 0),  c(0, -1, 0), c(0, 0, -1)),
  ci = c(a = -0.02, b = -5*10^(-3), d = -10^(-3), 
         a = -0.02, b = -5*10^(-3), d = -10^(-3)),
  method = "Nelder-Mead")

# Maximum likelihood estimate for effect of L1 density and L1 start
ML_3Pars_L1countL1width <- constrOptim(
  theta = c(a = ML_2Pars_L1count$par[1], ML_2Pars_L1count$par[2], 
            c = ML_L1width$par[2]),
  f = function(x) -AlleleFreqLogLik_4Par(
    Freqs = round(L1TotData$Freq[!blnNA], 0),
    Counts = rep(1, sum(!blnNA)),
    Predict = PredictMat[!blnNA, 1:3],
    a = x[1], b = x[2], c = x[3], d = 0, N = 10^4, 
    SampleSize = L1TotData$SampleSize[!blnNA],
    blnIns = L1TotData$blnIns[!blnNA], 
    LogRegCoeff = LogRegL1Ref$coefficients,
    DetectProb = L1TotData$DetectProb[!blnNA]),
  grad = NULL,
  ui = rbind(c(1, 0, 0),  c(0, 1, 0), c(0, 0, 1),     
             c(-1, 0, 0),  c(0, -1, 0), c(0, 0, -1)),
  ci = c(a = -0.02, b = -5*10^(-3), c = -10^(-6), 
         a = -0.02, b = -5*10^(-3), c = -10^(-6)),
  method = "Nelder-Mead")

# Maximum likelihood estimate for effect of L1 density, L1 start, and 
# full-length L1
ML_4Pars_L1countL1widthL1full <- constrOptim(
  theta = c(a = ML_2Pars_L1count$par[1], b = 0, 
            c = ML_L1widthL1full$par[2], d = ML_L1widthL1full$par[3]),
  f = function(x) -AlleleFreqLogLik_4Par(
    Freqs = round(L1TotData$Freq[!blnNA], 0),
    Counts = rep(1, sum(!blnNA)), 
    Predict = PredictMat[!blnNA, 1:3], 
    a = x[1], b = x[2], c = x[3], d = x[4], N = 10^4, 
    SampleSize = L1TotData$SampleSize[!blnNA],
    blnIns = L1TotData$blnIns[!blnNA], 
    LogRegCoeff = LogRegL1Ref$coefficients,
    DetectProb = L1TotData$DetectProb[!blnNA]),
  grad = NULL,
  ui = rbind(c(1, 0, 0, 0),  c(0, 1, 0, 0), c(0, 0, 1, 0),  c(0, 0, 0, 1),   
             c(-1, 0, 0, 0),  c(0, -1, 0, 0), c(0, 0, -1, 0),  c(0, 0, 0, -1)),
  ci = c(a = -0.01, b = -2*10^(-3), c = -10^(-6), d = -10^(-3), 
         a = -0.01, b = -2*10^(-3), c = -10^(-6), d = -10^(-3)),
  method = "Nelder-Mead")


###################################################
#                                                 #
#  Compare estimated and observed frequencies     #
#                                                 #
###################################################

# LogProbs <- AlleleFreqSampleProb(s = 0, N = 10^4, SampleSize = 2*2504)
# sum(is.infinite(LogProbs))
# length(LogProbs)
# min(LogProbs[!is.infinite(LogProbs)])
# idxFinite <- which(!is.infinite(LogProbs))
# plot(idxFinite, LogProbs[idxFinite])
# plot(idxFinite, LogProbs[idxFinite], xlim = c(4800, 5000))
# lchoose(5008, 1000)
#  k <- 200
#  SampleSize = 5008
#  integrate(function(x) AlleleFreqTime(x, s = 0, N = 10^4) * x^(k) *
#              (1 - x)^(SampleSize - k) , 0, 1)$value
#  integrate(function(x) log(AlleleFreqTime(x, s = 0, N = 10^4)) + 
#              k * log(x) + (SampleSize - k) * log(1 - x), 0, 1)$value
#  
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
Cols2Append <- t(sapply(list(ML_1Par, ML_L1width, ML_L1full, ML_2Pars_L1count, 
                             ML_L1widthL1full, ML_3Pars_L1countL1width,
                             ML_3Pars_L1countL1full
#         ML_4Pars_L1countL1widthL1full
         ), function(x){
           c(AIC = GetAIC(x), Pars = GetParVals(x))
         }))
# Combine AIC values into one vector
AICTab <- cbind(data.frame(
            NrParameters = c(1, 2, 2, 2, 3, 3, 3
#                             4
                             ),
            Predictor = c("none", "L1 start", "L1 full-length", "L1count",
                                  "L1 start and full-length",
                          "L1 count and L1 start",
                          "L1 count and L1 full"
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
                            a = x[1], b = x[2], c = 0, d = 0, N = 10^4,
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
                            a = x[1], b = 0, c = x[2], d = 0, N = 10^4,
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
                              a = x[1], b = 0, c = 0, d = x[2], N = 10^4,
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
                            a = x[1], b = x[2], c = x[3], d = 0, N = 10^4,
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
  a = x[1], b = x[2], c = x[3], d = x[4], N = 10^4,
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
                                    a = x[1], b = x[2], c = x[3], d = x[3], N = 10^4,
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
ProbL1 <- sapply(SVals, function(x) ProbAlleleIncluded(x,N = 10^4, SampleSize = 2*2504))
par(oma = c(7, 1, 0, 2), mfrow = c(2, 1), mai = c(0.5, 1, 0.5, 1))
plot(SCoeffVect, InsPerbp[2,], ylab = "LINE-1s per Mb", 
     xlab = "", ylim = c(0, 3), xlim = c(-0.0025, 0), main = "A")
text(SCoeffVect, InsPerbp[2,] + 2*10^(-1), names(SCoeffVect))
par(new = TRUE)
plot(SVals, ProbL1, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side = 4)
mtext("Inclusion probability", 4, line = 3)

# Plot expected frequency versus observed mean frequency
ExpL1 <- sapply(SVals, function(x) ExpAlleleFreq(x, N = 10^4, SampleSize = 2*2504))
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

# Create a vector of L1 start classes
L1TotData$L1widthClass <- cut(L1TotData$L1width, breaks = 
                                  seq(0, 7000, 1000))
MEInsCall$L1widthClass <- cut(MEInsCall$L1width, breaks = 
                                   seq(0, 7000, 1000))

MEInsCall$Freq
# Get mean L1 frequency per start
L1widthAggregated <- aggregate(L1TotData[,c("L1width", "L1Freq")], 
                               by = list(L1TotData$L1widthClass), 
                               FUN = function(x) mean(x, na.rm = T))
L1widthAggregated_Ins <- aggregate(MEInsCall[,c("L1width", "AF")], 
                               by = list(MEInsCall$L1widthClass), 
                               FUN = function(x) mean(x, na.rm = T))
plot(L1widthAggregated_Ins$L1width, L1widthAggregated_Ins$AF)

# Get sample size and create a range of s-values
SSize <- 2 * MEInsSamplesize
StartVals  <- seq(0, 6000, 100)
Full       <- StartVals == 6000
SVals <- ML_L1widthL1full$par[1] + ML_L1widthL1full$par[2]*StartVals +
  ML_L1widthL1full$par[3]*Full

# Plot expected frequency versus observed mean frequency
ExpL1width <- sapply(SVals, function(x) ExpAlleleFreq(x, N = 10^4, 
                                                      SampleSize = 2*MEInsSamplesize))
par( mfrow = c(1, 1))
plot(L1widthAggregated$L1width, 
     L1widthAggregated$L1Freq, xlab = "LINE-1 length",
     ylab = "Mean LINE-1 frequency")
lines(StartVals, ExpL1width )
mtext("Selection coefficient", 1, line = 3)
CreateDisplayPdf('D:/L1polymORF/Figures/FreqVsL1width_MELT.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7, width = 7)

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
                              a = x[1], b = 0, c = 0, d = 0, N = 10^4,
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
                                a = x[1], b = x[2], c = 0, d = 0, N = 10^4,
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
    a = x[1], b = 0, c = 0, d = 0, N = 10^4,
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
    a = x[1], b = x[2], c = 0, d = 0, N = 10^4,
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

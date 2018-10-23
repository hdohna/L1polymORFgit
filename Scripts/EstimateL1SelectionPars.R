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
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)

##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

# Specify file paths
DataPath        <- 'D:/L1polymORF/Data/'
G1000SamplePath <- 'D:/L1polymORF/Data/1000GenomeSampleInfo.txt'
GapPath         <- 'D:/L1polymORF/Data/Gap_hg19.txt'
L1GRPath        <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
L1RefRangePath  <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
ChrLPath        <- 'D:/L1polymORF/Data/ChromLengthsHg19.Rdata'
InputPath       <- 'D:/L1polymORF/Data/SingletonAnalysis_unphased.RData'
L1RefPath       <- 'D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv'
RegrOutputPath     <- "D:/L1polymORF/Data/L1RegressionResults.RData"
SelectTabOutPath    <- "D:/L1polymORF/Data/L1SelectionResults.csv"
SelectGenTabOutPath <- "D:/L1polymORF/Data/L1SelectionGeneResults.csv"
SelectResultOutPath <- "D:/L1polymORF/Data/L1SelectionResults.RData"

# Number of info columns in vcf file
NrInfoCols   <- 9

# False discovery rate for selected L1
FDR <- 0.1

# Specify range width for DNAse analysis
RangeWidth <- 10^6


##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("\n\nLoading and processing data ...")

# Create genomic range for MHC
MHC_GR <- GRanges(seqnames = "chr6",
                  IRanges(start = 28477797, end = 33448354))

# Load previously generated objects
load(InputPath)
load(L1GRPath)
load(ChrLPath)
load(L1RefRangePath)
load(RegrOutputPath)

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

# Table for each L1 how often it occurs in each population
UniquePops <- unique(SampleInfo$super_pop)
PopFreq <- sapply(UniquePops, function(x){
  blnPop <- Pops == x
  rowSums(L1_1000G[,SampleColumns[blnPop]])
})
colnames(PopFreq) <- UniquePops

# Match coefficients to 1000 genome data
ChromPosCoeff     <- paste(L1SingletonCoeffs$chromosome, L1SingletonCoeffs$Pos)
ChromPos1000G     <- paste(L1_1000G_reduced$chromosome, L1_1000G_reduced$POS)
MatchCoeff1000G   <- match(ChromPosCoeff, ChromPos1000G)
L1SingletonCoeffs <- cbind(L1SingletonCoeffs, PopFreq[MatchCoeff1000G,])

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

# Read in table with L1HS from the referemce genome
L1RefTab <- read.csv(L1RefPath)
L1RefGR <- makeGRangesFromDataFrame(L1RefTab, seqnames.field = "genoName",
                                    start.field = "genoStart",
                                    end.field = "genoEnd")
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
  L1SingletonCoeffs$L1End >= 6000
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
L1_1000G$blnOLGene  <- overlapsAny(L1_1000G_GR_hg19, GeneGR, ignore.strand = T)
L1_1000G$blnOLGeneSameStrand <- overlapsAny(L1_1000G_GR_hg19, GeneGR)
L1_1000G$blnOLProm     <- overlapsAny(L1_1000G_GR_hg19, PromGR, ignore.strand = T)
L1_1000G$blnOLExon     <- overlapsAny(L1_1000G_GR_hg19, ExonGR, ignore.strand = T)
L1_1000G$blnOLIntron   <- L1_1000G$blnOLGene & (!L1_1000G$blnOLExon)
L1_1000G$blnOLIntergen <- !(L1_1000G$blnOLGene | L1_1000G$blnOLProm)
L1_1000G$L1StartNum <- as.numeric(as.character(L1_1000G$L1Start))
L1_1000G$L1EndNum   <- as.numeric(as.character(L1_1000G$L1End))
L1_1000G$blnFull    <- L1_1000G$L1StartNum <= 1 & L1_1000G$L1EndNum >= 6000

# Create a variable indicating insertion type
L1_1000G$InsType <- "Intergenic"
L1_1000G$InsType[L1_1000G$blnOLProm]     <- "Promoter"
L1_1000G$InsType[L1_1000G$blnOLExon]     <- "Exon"
L1_1000G$InsType[L1_1000G$blnOLIntron]   <- "Intron"
table(L1_1000G$blnOLIntergen)

# Perform pairwise Wilcoxon test for differences in L1 frequencies
pairwise.wilcox.test(L1_1000G$Frequency, L1_1000G$InsType,
                     p.adjust.method = "BH")
# Average 
MeanFreqAgg <- aggregate(Frequency ~ InsType, data = L1_1000G, FUN = mean)
VarFreqAgg  <- aggregate(Frequency ~ InsType, data = L1_1000G, FUN = var)
L1_1000G$Dummy <- 1
NAgg  <- aggregate(Dummy ~ InsType, data = L1_1000G, FUN = sum)
StErr <- sqrt(VarFreqAgg$Frequency / NAgg$Dummy)

# Indicator variable for intersection with reference L1
blnOLGene_RefL1  <- overlapsAny(L1RefGR, GeneGR, ignore.strand = T)
blnOLGeneSameStrand_RefL1 <- overlapsAny(L1RefGR, GeneGR)
blnOLProm_RefL1   <- overlapsAny(L1RefGR, PromGR, ignore.strand = T)
blnOLExon_RefL1   <- overlapsAny(L1RefGR, ExonGR, ignore.strand = T)
blnOLIntron_RefL1 <- blnOLGene_RefL1 & (!blnOLExon_RefL1)

# Get number of insertions per bp
GeneTot     <- sum(width(GeneGR))
ExonTot     <- sum(width(ExonGR))
IntronTot   <- GeneTot - ExonTot
PromTot     <- sum(width(PromGR))
IntergenTot <- sum(as.numeric(ChromLengthsHg19)) - GeneTot - PromTot #- EnhancerTot

# Get mean frequency of L1 in different functional regions
MeanFreqs <- c(
               Promoter = mean(L1_1000G$Frequency[L1_1000G$blnOLProm]),
               Exon = mean(L1_1000G$Frequency[L1_1000G$blnOLExon]),
               Intron = mean(L1_1000G$Frequency[L1_1000G$blnOLIntron]),
               Intergenic = mean(L1_1000G$Frequency[L1_1000G$blnOLIntergen])
               )

# Plot distn of frequency of L1 in different functional regions
par(mfrow = c(1, 1))
hist(L1_1000G$Frequency[L1_1000G$blnOLProm], breaks = seq(0, 1, 0.01))
hist(L1_1000G$Frequency[L1_1000G$blnOLExon], breaks = seq(0, 1, 0.01))
hist(L1_1000G$Frequency[L1_1000G$blnOLIntron], breaks = seq(0, 1, 0.01))
hist(L1_1000G$Frequency[L1_1000G$blnOLIntergen], breaks = seq(0, 1, 0.01))
hist(sqrt(-log10(L1_1000G$Frequency[L1_1000G$blnOLProm])))
hist(-log10(L1_1000G$Frequency[L1_1000G$blnOLExon]))
hist(log10(L1_1000G$Frequency[L1_1000G$blnOLIntron]))
hist(log10(L1_1000G$Frequency[L1_1000G$blnOLIntergen]))

# Get number of L1 per Mb in different functional regions
InsPerbp <- 10^6 * rbind(
  c(
    Promoter = sum(blnOLProm_RefL1) / PromTot,
    Exon = sum(blnOLExon_RefL1) / ExonTot,
    Intron = sum(blnOLIntron_RefL1) / IntronTot,
    Intergenic = sum(!(blnOLGene_RefL1 | blnOLProm_RefL1)) / IntergenTot
    ),
  c(
    Promoter = sum(L1_1000G$blnOLProm) / PromTot,
    Exon = sum(L1_1000G$blnOLExon) / ExonTot,
    Intron = sum(L1_1000G$blnOLIntron) / IntronTot,
    Intergenic = sum(!(L1_1000G$blnOLGene | L1_1000G$blnOLProm)) / IntergenTot
    )
  
)
InsPerbp[1,] / InsPerbp[2,]
par(mfrow = c(2, 2), mai = c(1, 1, 1, 0))
bp1 <- barplot(InsPerbp, ylab = "L1 insertions per Mb", beside = T, 
        main = "A", las = 2, ylim = c(0, 2.5))
text(x = bp1[2,], y = InsPerbp[2,] + 0.2, # NAgg$Dummy[c(1, 5, 2, 4, 3)]
     NAgg$Dummy[c(4, 1, 3, 2)])
barplot(InsPerbp[2,] / InsPerbp[1,], ylab = "1000 genomes / reference",
        main = "B", las = 2)
bp1 <- barplot(MeanFreqs, ylab = "Mean frequency", main = "C", las = 2,
               ylim = c(0, 0.04))
AddErrorBars(MidX = bp1, MidY = MeanFreqs, ErrorRange = StErr,
             TipWidth = 0.1)
CreateDisplayPdf('D:/L1polymORF/Figures/PropL1InRegions.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7, width = 7)

###################################################
#                                                 #
#   Fit effect of insertion length on selection   #
#                                                 #
###################################################

# Match summary ranges to L1 ranges of 1000 genome data
L1SummaryOL <- findOverlaps(L1_1000G_GR_hg19, SummaryGR)
all(L1SummaryOL@from == 1:nrow(L1_1000G))
L1_1000G$L1Count <- DataPerSummaryGR$L1Count[L1SummaryOL@to]

# Get distance to nearest other L1
Dist2Nearest <- distanceToNearest(L1_1000G_GR_hg19)
L1_1000G$Dist2Nearest <-  Dist2Nearest@elementMetadata@listData$distance
max(L1_1000G$Dist2Nearest)

# Create a matrix of predictor variables (L1 start and boolean variable for)
PredictMat <- L1_1000G[, c("L1Count", "L1StartNum", "blnFull")]
blnNA <- sapply(1:nrow(L1_1000G), function(x) any(is.na(PredictMat[x,])))
sum(!blnNA)

# Plot log-likelihood for different selection coefficients
aVals <- seq(-0.002, 0.0003, 0.0001)
LikVals <- sapply(aVals, function(x) {
  print(x)
  LL_FPrime = AlleleFreqLogLik_4Par(
    Freqs = (L1_1000G$Frequency * 2*2504)[!blnNA],
    Counts = rep(1, sum(!blnNA)),
    Predict = PredictMat[!blnNA,],
    a = x, b = 0, c = 0, d = 0, N = 10^4,
    SampleSize = 2*2504, blnUseFPrime = T)
  })
par(mfrow = c(1, 1))
plot(aVals, LikVals, type = "l", col = "red")

# Estimate maximum likelihood for a single selection coefficient
ML_1Par <-  constrOptim(theta = c(a = -0.0004),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = (L1_1000G$Frequency * 2*2504)[!blnNA],
                            Counts = rep(1, sum(!blnNA)),
                            Predict = PredictMat[!blnNA,],
                            a = x[1], b = 0, c = 0, d = 0, N = 10^4,
                            SampleSize = 2*2504),
                          grad = NULL,
                          ui = rbind(1,-1),
                          ci = c(a = -0.01, a = -0.01),
                          method = "Nelder-Mead")
cat("done!\n")

# Get maximum likelihood estimate for effect of L1 start on selection
cat("Estimate effect of L1 start on selections ...")
ML_L1start <-  constrOptim(theta = c(a = ML_1Par$par, c = 0),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = (L1_1000G$Frequency * 2*2504)[!blnNA], 
                            Counts = rep(1, sum(!blnNA)), 
                            Predict = PredictMat[!blnNA,], 
                            a = x[1], b = 0, c = x[2], d = 0, N = 10^4, 
                            SampleSize = 2*2504),
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
                      Freqs = (L1_1000G$Frequency * 2*2504)[!blnNA], 
                      Counts = rep(1, sum(!blnNA)), 
                      Predict = PredictMat[!blnNA,], 
                      a = x[1], b = 0, c = 0, d = x[2], N = 10^4, 
                      SampleSize = 2*2504),
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
ML_L1startL1full <- constrOptim(theta = c(a = ML_L1start$par[1], b = ML_L1start$par[2], 
                                          c = ML_L1full$par[2]),
                  f = function(x) -AlleleFreqLogLik_4Par(
                    Freqs = (L1_1000G$Frequency * 2*2504)[!blnNA], 
                    Counts = rep(1, sum(!blnNA)), 
                    Predict = PredictMat[!blnNA,], 
                    a = x[1], b = 0, c = x[2], d = x[3], N = 10^4, 
                    SampleSize = 2*2504),
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
                            Freqs = (L1_1000G$Frequency * 2*2504)[!blnNA], 
                            Counts = rep(1, sum(!blnNA)), 
                            Predict = PredictMat[!blnNA,], 
                            a = x[1], b = x[2], c = 0, d = 0, N = 10^4, 
                            SampleSize = 2*2504),
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
    Freqs = (L1_1000G$Frequency * 2*2504)[!blnNA], 
    Counts = rep(1, sum(!blnNA)), 
    Predict = PredictMat[!blnNA,], 
    a = x[1], b = x[2], c = 0, d = x[3], N = 10^4, 
    SampleSize = 2*2504),
  grad = NULL,
  ui = rbind(c(1, 0, 0),  c(0, 1, 0), c(0, 0, 1),     
             c(-1, 0, 0),  c(0, -1, 0), c(0, 0, -1)),
  ci = c(a = -0.02, b = -5*10^(-3), d = -10^(-3), 
         a = -0.02, b = -5*10^(-3), d = -10^(-3)),
  method = "Nelder-Mead")

# Maximum likelihood estimate for effect of L1 density and L1 start
ML_3Pars_L1countL1start <- constrOptim(
  theta = c(a = ML_2Pars_L1count$par[1], ML_2Pars_L1count$par[2], 
            c = ML_L1start$par[2]),
  f = function(x) -AlleleFreqLogLik_4Par(
    Freqs = (L1_1000G$Frequency * 2*2504)[!blnNA], 
    Counts = rep(1, sum(!blnNA)), 
    Predict = PredictMat[!blnNA,], 
    a = x[1], b = x[2], c = x[3], d = 0, N = 10^4, 
    SampleSize = 2*2504),
  grad = NULL,
  ui = rbind(c(1, 0, 0),  c(0, 1, 0), c(0, 0, 1),     
             c(-1, 0, 0),  c(0, -1, 0), c(0, 0, -1)),
  ci = c(a = -0.02, b = -5*10^(-3), c = -10^(-6), 
         a = -0.02, b = -5*10^(-3), c = -10^(-6)),
  method = "Nelder-Mead")

# Maximum likelihood estimate for effect of L1 density, L1 start, and 
# full-length L1
ML_4Pars_L1countL1startL1full <- constrOptim(
  theta = c(a = ML_2Pars_L1count$par[1], b = ML_2Pars_L1count$par[2], 
            c = ML_L1start$par[2], d = ML_L1full$par[2]),
  f = function(x) -AlleleFreqLogLik_4Par(
    Freqs = (L1_1000G$Frequency * 2*2504)[!blnNA], 
    Counts = rep(1, sum(!blnNA)), 
    Predict = PredictMat[!blnNA,], 
    a = x[1], b = x[2], c = x[3], d = x[4], N = 10^4, 
    SampleSize = 2*2504),
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
Cols2Append <- t(sapply(list(ML_1Par, ML_L1start, ML_L1full, ML_2Pars_L1count, 
                             ML_L1startL1full, ML_3Pars_L1countL1start,
                             ML_3Pars_L1countL1full,
         ML_4Pars_L1countL1startL1full), function(x){
           c(AIC = GetAIC(x), Pars = GetParVals(x))
         }))
# Combine AIC values into one vector
AICTab <- cbind(data.frame(
            NrParameters = c(1, 2, 2, 2, 3, 3, 3, 4),
            Predictor = c("none", "L1 start", "L1 full-length", "L1count",
                                  "L1 start and full-length",
                          "L1 count and L1 start",
                          "L1 count and L1 full",
                          "L1 start, L1 full-length, L1count"),
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
PredictMatGeneOL <- L1_1000G[, c("blnOLExon", "blnOLIntron", "blnOLProm")]
PredictMatGeneOL2 <- L1_1000G[, c("blnOLGene", "blnOLIntron", "blnOLProm")]

# Estimate maximum likelihood for a single selection coefficient
ML_1Par_all <-  constrOptim(theta = c(a = -0.0004),
                        f = function(x) -AlleleFreqLogLik_4Par(
                          Freqs = (L1_1000G$Frequency * 2*2504),
                          Counts = rep(1, nrow(PredictMatGeneOL)),
                          Predict = PredictMatGeneOL,
                          a = x[1], b = 0, c = 0, d = 0, N = 10^4,
                          SampleSize = 2*2504),
                        grad = NULL,
                        ui = rbind(1,-1),
                        ci = c(a = -0.01, a = -0.01),
                        method = "Nelder-Mead")

# Get maximum likelihood estimate for effect of exonic L1 on selection
cat("Estimate effect of exon overlap on selections ...")
ML_L1Exon <-  constrOptim(theta = c(a = ML_1Par_all$par, b = 0),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = (L1_1000G$Frequency * 2*2504),
                            Counts = rep(1, nrow(PredictMatGeneOL)),
                            Predict = PredictMatGeneOL,
                            a = x[1], b = x[2], c = 0, d = 0, N = 10^4,
                            SampleSize = 2*2504),
                          grad = NULL,
                           ui = rbind(c(1, 0),  c(0, 1),   
                                      c(-1, 0), c(0, -1)),
                           ci = c(a = -0.01, b = -10^(-2), 
                                  a = -0.01, b = -10^(-2)),
                           method = "Nelder-Mead")
cat("done!\n")

# Get maximum likelihood estimate for effect of exonic L1 on selection
cat("Estimate effect of exon or intron overlap on selections ...")
ML_L1ExonOrIntron <-  constrOptim(theta = c(a = ML_1Par_all$par, b = 0),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = (L1_1000G$Frequency * 2*2504),
                            Counts = rep(1, nrow(PredictMatGeneOL2)),
                            Predict = PredictMatGeneOL,
                            a = x[1], b = x[2], c = 0, d = 0, N = 10^4,
                            SampleSize = 2*2504),
                          grad = NULL,
                          ui = rbind(c(1, 0),  c(0, 1),   
                                     c(-1, 0), c(0, -1)),
                          ci = c(a = -0.01, b = -10^(-2), 
                                 a = -0.01, b = -10^(-2)),
                          method = "Nelder-Mead")
cat("done!\n")

# Get maximum likelihood estimate for effect of intronic L1 on selection
cat("Estimate effect of intron overlap on selections ...")
ML_L1Intron <-  constrOptim(theta = c(a = ML_1Par_all$par, c = 0),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = (L1_1000G$Frequency * 2*2504),
                            Counts = rep(1, nrow(PredictMatGeneOL)),
                            Predict = PredictMatGeneOL,
                            a = x[1], b = 0, c = x[2], d = 0, N = 10^4,
                            SampleSize = 2*2504),
                          grad = NULL,
                          ui = rbind(c(1, 0),  c(0, 1),   
                                     c(-1, 0), c(0, -1)),
                          ci = c(a = -0.01, c = -10^(-2), 
                                 a = -0.01, c = -10^(-2)),
                          method = "Nelder-Mead")
cat("done!\n")

# Get maximum likelihood estimate for effect of intronic L1 on selection
cat("Estimate effect of promoter overlap on selections ...")
ML_L1Prom <-  constrOptim(theta = c(a = ML_1Par_all$par, d = 0),
                            f = function(x) -AlleleFreqLogLik_4Par(
                              Freqs = (L1_1000G$Frequency * 2*2504),
                              Counts = rep(1, nrow(PredictMatGeneOL)),
                              Predict = PredictMatGeneOL,
                              a = x[1], b = 0, c = 0, d = x[2], N = 10^4,
                              SampleSize = 2*2504),
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
                          theta = c(a = ML_1Par_all$par, 
                                    b = ML_L1Exon$par[2],
                                    c = ML_L1Intron$par[2]),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = (L1_1000G$Frequency * 2*2504),
                            Counts = rep(1, nrow(PredictMatGeneOL)),
                            Predict = PredictMatGeneOL,
                            a = x[1], b = x[2], c = x[3], d = 0, N = 10^4,
                            SampleSize = 2*2504),
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
    Freqs = (L1_1000G$Frequency * 2*2504),
    Counts = rep(1, nrow(PredictMatGeneOL)),
    Predict = PredictMatGeneOL,
    a = x[1], b = x[2], c = x[3], d = 0, N = 10^4,
    SampleSize = 2*2504),
  grad = NULL,
  ui = rbind(c(1, 0, 0, 0),  c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1),  
             c(-1, 0, 0, 0), c(0, -1, 0, 0) , c(0, 0, -1, 0), c(0, 0, 0, -1)),
  ci = c(a = -0.01, b = -10^(-2), c = -10^(-2), d = -10^(-2), 
         a = -0.01, b = -10^(-2), c = -10^(-2), d = -10^(-2)),
  method = "Nelder-Mead")
cat("done!\n")

# Get columns of AIC and parameter values
Cols2Append <- t(sapply(list(ML_1Par_all, ML_L1Exon, ML_L1Intron, ML_L1Prom, 
                             ML_L1ExonOrIntron, ML_L1ExonIntron,
                             ML_L1ExonIntronProm), 
                        function(x){
                               c(NrParameters = GetNPar(x), AIC = GetAIC(x), 
                                 Pars = GetParVals(x))
                             }))
# Combine AIC values into one vector
AICTabGene <- cbind(data.frame(
  Predictor = c("none", "Exon", "Intron", "Promoter",
                "Exon or intron", "Exon and intron", 
                "Exon, intron, and promoter"),
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
SSize <- 2*2504
SVals  <- seq(-0.0025, -0.00001, 0.00001)

# Plot probability for inclusion versus number of LINE-1 per Mb
ProbL1 <- sapply(SVals, function(x) ProbAlleleIncluded(x,N = 10^4, SampleSize = 2*2504))
par(oma = c(7, 1, 0, 2), mfrow = c(2, 1), mai = c(0.5, 1, 0.5, 1))
plot(SCoeffVect, InsPerbp[2,], ylab = "LINE-1s per Mb", 
     xlab = "", ylim = c(0, 2.5), xlim = c(-0.0025, 0), main = "A")
text(SCoeffVect, InsPerbp[2,] + 2*10^(-1), names(SCoeffVect))
par(new = TRUE)
plot(SVals, ProbL1, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side = 4)
mtext("Inclusion probability", 4, line = 3)

# Plot expected frequency versus observed mean frequency
ExpL1 <- sapply(SVals, function(x) ExpAlleleFreq(x, N = 10^4, SampleSize = 2*2504))
plot(SCoeffVect, MeanFreqs*SSize, ylab = "Mean LINE-1 frequency", 
     xlab = "", xlim = c(-0.0025, 0.0001), ylim = c(20, 155), main = "B")
text(SCoeffVect + c(0.0002, 0, -0.0001, -0.0002), MeanFreqs*SSize + 10, names(SCoeffVect))
lines(SVals, ExpL1)
mtext("Selection coefficient", 1, line = 3)
CreateDisplayPdf('D:/L1polymORF/Figures/SelectionPerRegion.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7, width = 7)

###################################################
#                                                 #
#   Plot frequency vs. insertion length           #
#                                                 #
###################################################

# Create a vector of L1 start classes
L1_1000G$L1StartNumClass <- cut(L1_1000G$L1StartNum, breaks = 
                                  seq(0, 6100, 500))
L1_1000G$Frequency

# Get mean L1 frequency per start
L1StartAggregated <- aggregate(L1_1000G[,c("L1StartNum", "Frequency")], 
                               by = list(L1_1000G$L1StartNumClass), FUN = mean)

# Get sample size and create a range of s-values
SSize <- 2*2504
StartVals  <- seq(0, 6000, 100)
Full       <- StartVals == 0
SVals <- ML_L1startL1full$par[1] + ML_L1startL1full$par[2]*StartVals +
  ML_L1startL1full$par[3]*Full

# Plot expected frequency versus observed mean frequency
ExpL1Start <- sapply(SVals, function(x) ExpAlleleFreq(x, N = 10^4, SampleSize = 2*2504))
par( mfrow = c(1, 1))
plot(L1StartAggregated$L1StartNum, 
     L1StartAggregated$Frequency * SSize, xlab = "5' start of LINE-1",
     ylab = "Mean LINE-1 frequency")
lines(StartVals, ExpL1Start)
mtext("Selection coefficient", 1, line = 3)
CreateDisplayPdf('D:/L1polymORF/Figures/FreqVsL1Start.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7, width = 7)

###################################################
#                                                 #
#   Fit effect of strandedness on selection       #
#                                                 #
###################################################

# Create a matrix of predictor variables (L1 start and boolean variable for)
PredictMatWithinGene <- L1_1000G[L1_1000G$blnOLGene, 
                                 c( "blnOLGeneSameStrand", "blnOLGene", "blnOLGene")]

# Estimate maximum likelihood for a single selection coefficient
ML_1Par_gene <- constrOptim(theta = c(a = -0.0004),
                            f = function(x) -AlleleFreqLogLik_4Par(
                              Freqs = (L1_1000G$Frequency[L1_1000G$blnOLGene] * 2*2504),
                              Counts = rep(1, sum(L1_1000G$blnOLGene)),
                              Predict = PredictMatWithinGene,
                              a = x[1], b = 0, c = 0, d = 0, N = 10^4,
                              SampleSize = 2*2504),
                            grad = NULL,
                            ui = rbind(1,-1),
                            ci = c(a = -0.01, a = -0.01),
                            method = "Nelder-Mead")

# Get maximum likelihood estimate for effect of exonic L1 on selection
cat("Estimate effect of same strand overlap on selections ...")
ML_L1SameStrand <-  constrOptim(theta = c(a = ML_1Par_gene$par, b = 0),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = (L1_1000G$Frequency[L1_1000G$blnOLGene] * 2*2504),
                            Counts = rep(1, sum(L1_1000G$blnOLGene)),
                            Predict = PredictMatWithinGene,
                            a = x[1], b = x[2], c = 0, d = 0, N = 10^4,
                            SampleSize = 2*2504),
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


###################################################
#                                                 #
#   Fit effect of singleton coef. on selection    #
#                                                 #
###################################################

# Create a matrix of predictor variables 
PredictMat <- L1SingletonCoeffs[, c("coef", "L1Start", "L1Start")]
blnNA <- sapply(1:nrow(L1SingletonCoeffs), function(x) any(is.na(PredictMat[x,])))

# Determine maximum likelihood with one parameter (selection coefficient)
cat("Maximizing likelihood for one parameter (selection coefficient) ...")
ML_1Par_coef <- optim(par = c(a = ML_1Par$par),
                      fn = function(x) -AlleleFreqLogLik_4Par(
                        Freqs = (L1SingletonCoeffs$Freq * 2*2504)[!blnNA], 
                        Counts = rep(1, sum(!blnNA)), 
                        Predict = PredictMat[!blnNA,], 
                        a = x, b = 0, c = 0, d = 0, N = 10^4, 
                        SampleSize = 2*2504),
                      lower = -0.01, upper = 0.01,
                      method = "L-BFGS-B")
ML_2Pars_L1coef <- constrOptim(
  theta = c(a = ML_1Par_coef$par, b = 0),
  f = function(x) -AlleleFreqLogLik_4Par(
    Freqs = (L1SingletonCoeffs$Freq * 2*2504)[!blnNA], 
    Counts = rep(1, sum(!blnNA)), 
    Predict = PredictMat[!blnNA,], 
    a = x[1], b = x[2], c = 0, d = 0, N = 10^4, 
    SampleSize = 2*2504),
  grad = NULL,
  ui = rbind(c(1, 0),  c(0, 1),  
             c(-1, 0), c(0, -1)),
  ci = c(a = -0.01, b = -2*10^(-3), 
         a = -0.01, b = -2*10^(-3)),
  method = "Nelder-Mead")
cat("done!\n")

# Save everything
save.image(SelectResultOutPath)

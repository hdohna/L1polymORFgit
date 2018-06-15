# The script below analyzes overlap of L1 with L1
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
library(KernSmooth)
library(glmnet)
library(org.Hs.eg.db)

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

# Number of info columns in vcf file
NrInfoCols   <- 9


##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("Loading and processing data ...")

# Load previously generated objects
load(L1GRPath)
load(ChrLPath)
load(L1RefRangePath)
load(InputPath)

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
PromGR <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene)
cat("done!\n")

##########################################
#                                        #
#        Add columns                     #
#                                        #
##########################################

# Turn factors into numeric values
L1SingletonCoeffs$L1Start <- as.numeric(as.character(L1SingletonCoeffs$L1Start))
L1SingletonCoeffs$L1End <- as.numeric(as.character(L1SingletonCoeffs$L1End))

# Indicator for full-length
L1SingletonCoeffs$blnFull <- L1SingletonCoeffs$L1Start <= 3 &
  L1SingletonCoeffs$L1End >= 6000
sum(L1SingletonCoeffs$InsLength <= 100)

# Indicator for significant effect
L1SingletonCoeffs$blnSig <- p.adjust(L1SingletonCoeffs$Pr...z..) < 0.05
hist(L1SingletonCoeffs$Pr...z.., breaks = seq(0, 1, 0.005))

# Indicator for positive selection
L1SingletonCoeffs$blnSelect <- L1SingletonCoeffs$blnSig &
  L1SingletonCoeffs$coef < 0

# Indicator for negative selection
L1SingletonCoeffs$blnNotSelect <- L1SingletonCoeffs$blnSig &
  L1SingletonCoeffs$coef > 0
sum(L1SingletonCoeffs$blnNotSelect)
sum(L1SingletonCoeffs$blnSelect)

# Indicator fo selection (+1 = positive, -1 = negative, 0 = neutral)
L1SingletonCoeffs$SelectInd <- 0
L1SingletonCoeffs$SelectInd[L1SingletonCoeffs$blnSelect]     <- 1
L1SingletonCoeffs$SelectInd[L1SingletonCoeffs$blnNotSelect]  <- -1

# Caclulate distance to genes
L1SingletonCoeffs$Dist2Gene <- Dist2Closest(L1SingletonCoeffs_GR, GRgenes_hg19)

# Caclulate logarithm of distance to genes
L1SingletonCoeffs$LogDist2Gene <- log(L1SingletonCoeffs$Dist2Gene + 0.1)

##########################################
#                                        #
#     Regress coefficients               #
#                                        #
##########################################

# Form a subset of coefficients with nonzero standard error
L1SinglCoeff_nonzeroSE <- subset(L1SingletonCoeffs, subset = se.coef. > 0)

# Regress coefficients vs frequency
par(mfrow = c(1, 1))
plot(L1SingletonCoeffs$Freq, L1SingletonCoeffs$coef, 
     xlab = "Insertion length", ylab = "Singleton coefficient",
     xlim = c(0, 0.01))
L1CoefVsFreqSmoothed <- supsmu(L1SinglCoeff_nonzeroSE$Freq, L1SinglCoeff_nonzeroSE$coef,
                               wt = 1/L1SinglCoeff_nonzeroSE$se.coef.)
lines(L1CoefVsFreqSmoothed$x, L1CoefVsFreqSmoothed$y, col = "red")
lines(c(0, 1), c(0, 0), col = "blue")
LML1CoefVsFreq <- lm(coef ~ Freq, data = L1SinglCoeff_nonzeroSE, weights = 1/se.coef.)
summary(LML1CoefVsFreq)

# Regress coefficients vs L1 start
par(mfrow = c(1, 1))
plot(L1SingletonCoeffs$L1Start, L1SingletonCoeffs$coef, 
     xlab = "L1 start", ylab = "Singleton coefficient")
L1CoefVsL1StartSmoothed <- supsmu(L1SinglCoeff_nonzeroSE$L1Start, 
                                  L1SinglCoeff_nonzeroSE$coef,
                               wt = 1/L1SinglCoeff_nonzeroSE$se.coef.)
lines(L1CoefVsL1StartSmoothed$x, L1CoefVsL1StartSmoothed$y, col = "red")
lines(c(0, 10^4), c(0, 0), col = "blue")
LML1CoefVsL1Start <- lm(coef ~ L1Start, data = L1SinglCoeff_nonzeroSE, weights = 1/se.coef.)
summary(LML1CoefVsL1Start)

##########################################
#                                        #
#   Regress intersection with genes      #
#                                        #
##########################################

#######
# Regress against L1 start
#######

# Indicator variable for intersection with genes
L1_1000G$blnOLGene  <- overlapsAny(L1_1000G_GR_hg19, GRgenes_hg19)
L1_1000G$blnOLProm  <- overlapsAny(L1_1000G_GR_hg19, PromGR)
L1_1000G$L1StartNum <- as.numeric(as.character(L1_1000G$L1Start))
L1_1000G$L1EndNum   <- as.numeric(as.character(L1_1000G$L1End))
L1_1000G$blnFull    <- L1_1000G$L1StartNum <= 1 & L1_1000G$L1EndNum >= 6000
hist(L1_1000G$L1StartNum, breaks = 0:6100)
sum(L1_1000G$L1StartNum == 2, na.rm = T)
sum(L1_1000G$blnOLProm)

# Create transparent point color
PCol <- rgb(0,0,0, alpha = 0.1)

# Regress indicator for gene overlap against L1 start abd frequency
L1_1000G_subset <- L1_1000G[L1_1000G$Frequency >= 0.01, ]
L1_1000G_subset <- L1_1000G
LogRegGeneOLVsL1Start <- glm(blnOLGene ~ L1StartNum + Frequency, 
                             family = "binomial", data = L1_1000G_subset)
summary(LogRegGeneOLVsL1Start)

# Regress indicator for promoter overlap against L1 start abd frequency
LogRegPromOLVsL1Start <- glm(blnOLProm ~ L1StartNum  + Frequency, 
                             family = "binomial", data = L1_1000G)
summary(LogRegPromOLVsL1Start)

# Plot proportion gene overlap against L1 start
L1OLVsL1StartSmoothed <- supsmu(L1_1000G_subset$L1StartNum,  1*L1_1000G_subset$blnOLGene)
plot(L1OLVsL1StartSmoothed$x, L1OLVsL1StartSmoothed$y, type = "l", xlab = "L1 start",
     ylab = "Proportion of L1s in genes")
points(L1_1000G$L1StartNum[blnOLGene], rep(0.31, sum(blnOLGene)), col = PCol, 
       pch = 16)
plot(L1OLVsL1StartSmoothed$x, L1OLVsL1StartSmoothed$y, type = "l", xlab = "L1 start",
     ylab = "Proportion of L1s in genes", xlim = c(0, 40),
     ylim = c(0.38, 0.55))
points(L1_1000G$L1StartNum[blnOLGene], rep(0.4, sum(blnOLGene)), col = PCol, 
       pch = 16)

# Plot proportion gene overlap against L1 frequency
L1OLVsFreqSmoothed <- supsmu(L1_1000G$Frequency,  1*blnOLGene)
plot(L1OLVsFreqSmoothed$x, L1OLVsFreqSmoothed$y, type = "l")
points(L1_1000G$Frequency[blnOLGene], rep(0.4, sum(blnOLGene)))
plot(L1OLVsFreqSmoothed$x, L1OLVsFreqSmoothed$y, type = "l",
     xlim = c(0, 0.01))
points(L1_1000G$Frequency[blnOLGene], rep(0.4, sum(blnOLGene)), col = PCol, 
       pch = 16)

# Frequency and start range with high proportion in genes
blnLowFreq  <- L1_1000G$Frequency >= 0.001 & L1_1000G$Frequency <= 0.003
blnLowStart <- L1_1000G$L1StartNum >= 10 & L1_1000G$L1StartNum <= 15
sum(blnLowFreq & blnLowStart, na.rm = T)
mean(L1_1000G$blnOLGene[blnLowFreq & blnLowStart], na.rm = T)
mean(L1_1000G$blnOLGene[blnLowFreq], na.rm = T)
mean(L1_1000G$blnOLGene[blnLowStart], na.rm = T)

##########################################
#                                        #
#   Export IDs of intersecting genes     #
#                                        #
##########################################

# Define a cut-off value for L1 allele frequency
FreqCutOff <- 0.005
blnLowFreq <- L1_1000G$Frequency <= FreqCutOff

# Get genes that overlap with low and high frequency L1a
GeneGR_All      <- subsetByOverlaps(GRgenes_hg19, L1_1000G_GR_hg19)
GeneGR_LowFreq  <- subsetByOverlaps(GRgenes_hg19, L1_1000G_GR_hg19[blnLowFreq])
GeneGR_HiFreq   <- subsetByOverlaps(GRgenes_hg19, L1_1000G_GR_hg19[!blnLowFreq])

# Check whether gene size determines overlap
blnGeneOL <- overlapsAny(GRgenes_hg19, L1_1000G_GR_hg19)
t.test(log(width(GRgenes_hg19)) ~ blnGeneOL)
boxplot(log(width(GRgenes_hg19)) ~ blnGeneOL)

# Write out gene ids for different genelists
writeLines(GeneGR_All@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_All")
writeLines(GeneGR_LowFreq@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_LowFreq")
writeLines(GeneGR_HiFreq@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_HiFreq")

##########################################
#                                                 #
#   Sample genes at random according to width     #
#                                                 #
##########################################

# Sampled gene indices
idxSampledGenes <- sample(length(GRgenes_hg19), 
                         sum(L1_1000G$blnOLGene, na.rm = T), prob = width(GRgenes_hg19))
SampledGeneGR <- GRgenes_hg19[idxSampledGenes]
writeLines(SampledGeneGR@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/SampledGeneIDs")

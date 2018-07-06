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
library(UniProt.ws)
library(gee)
library(Homo.sapiens)
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
#UniProtPath     <- 'D:/L1polymORF/Data/UniProtKeywordsGO.txt'
UniProtPath     <- 'D:/L1polymORF/Data/UniProtKeywordsInterPro.txt'
InterProPath    <- 'D:/L1polymORF/Data/InterPro.txt'
DNAseDataPath   <- 'D:/L1polymORF/Data/DNAseInfo.RData'
GeneExpPath     <- 'D:/L1polymORF/Data/gtexGene.txt'
GExpTissuePath  <- 'D:/L1polymORF/Data/gtexTissue.txt'
InputPath       <- 'D:/L1polymORF/Data/SingletonAnalysis_unphased.RData'
L1RefPath       <- 'D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv'
CpGPath         <- 'D:/L1polymORF/Data/CpG_hg19.txt'
ParalogPath     <- "D:/L1polymORF/Data/Paralogs_hg19.RData"

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

cat("Loading and processing data ...")

# Create genomic range for MHC
MHC_GR <- GRanges(seqnames = "chr6",
                  IRanges(start = 28477797, end = 33448354))

# Load previously generated objects
load(L1GRPath)
load(ChrLPath)
load(L1RefRangePath)
load(InputPath)
load(ParalogPath)

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

# Read UniProt data
UniProtData  <- read.delim(UniProtPath)
UniProtData  <- UniProtData[,1:5]
InterProData <- read.delim(InterProPath)

# Read DNAse hypersensitivity data (generated in script CombineEncodeInfo_DNAse)
load(DNAseDataPath)

# Get gene expression data
GExpData <- read.delim(GeneExpPath, header = F, 
              col.names = c("chrom", "start", "end", "name",
                  "score", "strand", "geneId", "geneType", "expCount", 
                  "expScores"))
GExpGR <- makeGRangesFromDataFrame(GExpData)
GExpByTissue <- t(sapply(as.character(GExpData$expScores), function(x){
  strsplit(x, ",")[[1]]
}))
GExpByTissue <- as.data.frame(GExpByTissue)
for (i in 1:ncol(GExpByTissue)){
  GExpByTissue[[i]] <- as.numeric(GExpByTissue[[i]])
}
GexpTissue   <- read.delim(GExpTissuePath, header = F, 
     col.names = c("id", "name", "description", "organ", "color"))
colnames(GExpByTissue) <- GexpTissue$name

# Read and process table with regulatory elements
# RegTable <- read.table("D:/L1polymORF/Data/EncodeBroadHMM/ChromHMMcombined.txt",
#                             header = T)
# idxHetero   <- grep("Heterochrom", RegTable$name)
# RegGR       <- makeGRangesFromDataFrame(RegTable)
# HeteroGR    <- RegGR[idxHetero]

# Read in table with L1HS from the referemce genome
L1RefTab <- read.csv(L1RefPath)
L1RefGR <- makeGRangesFromDataFrame(L1RefTab, seqnames.field = "genoName",
                                    start.field = "genoStart",
                                    end.field = "genoEnd")

# Read in CpG data and turn it into GRanges
CpGtable <- read.delim(CpGPath, header = T)
CpGGR    <- makeGRangesFromDataFrame(CpGtable, start.field = "chromStart",
                                     end.field = "chromEnd")
cat("done!\n")

##########################################
#                                        #
#        Add columns                     #
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
L1SingletonCoeffs$blnNotSelect <- L1SingletonCoeffs$blnSig &
  L1SingletonCoeffs$coef > 0
sum(L1SingletonCoeffs$blnNotSelect)
sum(L1SingletonCoeffs$blnSelect)

# Indicator fo selection (+1 = positive, -1 = negative, 0 = neutral)
L1SingletonCoeffs$SelectInd <- 0
L1SingletonCoeffs$SelectInd[L1SingletonCoeffs$blnSelect]     <- 1
L1SingletonCoeffs$SelectInd[L1SingletonCoeffs$blnNotSelect]  <- -1

# Caclulate distance to genes
L1SingletonCoeffs$Dist2Gene <- Dist2Closest(L1SingletonCoeffs_GR, GeneGR)
L1SingletonCoeffs$blnOLGene <- L1SingletonCoeffs$Dist2Gene == 0

# Caclulate logarithm of distance to genes
L1SingletonCoeffs$LogDist2Gene <- log(L1SingletonCoeffs$Dist2Gene + 0.1)

# Standardize SE ratio to one
L1SingletonCoeffs$SE_RatioSt <- 1/L1SingletonCoeffs$se.coef./
  mean(1/L1SingletonCoeffs$se.coef.)


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

##########################################
#                                        #
#     Regress coefficients               #
#                                        #
##########################################

cat("Regress coefficients ...")

# Form a subset of coefficients with nonzero standard error
L1SinglCoeff_nonzeroSE <- subset(L1SingletonCoeffs, subset = se.coef. > 0)

# Regress coefficients vs frequency
par(mfrow = c(1, 1))
plot(L1SingletonCoeffs$Freq, L1SingletonCoeffs$coef, 
     xlab = "Frequency", ylab = "Singleton coefficient",
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

# Regress p-values vs L1 start
par(mfrow = c(1, 1))
plot(L1SingletonCoeffs$L1Start, log10(L1SingletonCoeffs$Pr...z..), 
     xlab = "L1 start", ylab = "Singleton coefficient")
L1PVsL1StartSmoothed <- supsmu(L1SinglCoeff_nonzeroSE$L1Start, 
                                  log10(L1SinglCoeff_nonzeroSE$Pr...z..),
                                  wt = 1/L1SinglCoeff_nonzeroSE$se.coef.)
lines(L1PVsL1StartSmoothed$x, L1PVsL1StartSmoothed$y, col = "red")
LML1PVsL1Start <- glm(blnSig ~ L1Start + Freq, data = L1SinglCoeff_nonzeroSE, 
                      weights = 1/se.coef., family = "quasibinomial")
summary(LML1PVsL1Start)
LML1PVsFreq <- glm(blnSig ~ Freq, data = L1SinglCoeff_nonzeroSE, 
                      weights = 1/se.coef., family = "quasibinomial")
summary(LML1PVsFreq)

# Check whether coefficients differ by L1 inside or outside genes
t.test(coef ~ blnOLGene, data = L1SingletonCoeffs)

cat("done!\n")

##########################################
#                                        #
#     Regress selection indicator        #
#                                        #
##########################################

cat("Regressing selection indicator ...")

# Determine whether indicator of negative selection depends on insertion length
# or frequency

LML1NegSelVsL1StartFreq <- glm(blnNotSelect ~ Freq + L1Start, 
                           data = L1SinglCoeff_nonzeroSE, 
                           weights = 1/se.coef., family = "quasibinomial")
summary(LML1NegSelVsL1StartFreq)
LML1NegSelVsL1Start <- glm(blnNotSelect ~ L1Start, 
                           data = L1SinglCoeff_nonzeroSE, 
                      weights = 1/se.coef., family = "quasibinomial")
summary(LML1NegSelVsL1Start)

LML1NegSelVsFreq <- glm(blnNotSelect ~ Freq, data = L1SinglCoeff_nonzeroSE, 
                   weights = 1/se.coef., family = "quasibinomial")
summary(LML1NegSelVsFreq)

# Determine whether L1s with negative selection signal are more likely
# to be in vs out of genes
table(L1SingletonCoeffs$Dist2Gene == 0, L1SingletonCoeffs$blnNotSelect)
fisher.test(L1SingletonCoeffs$Dist2Gene == 0, L1SingletonCoeffs$blnNotSelect)
chisq.test(L1SingletonCoeffs$Dist2Gene == 0, L1SingletonCoeffs$blnNotSelect)

# Determine whether L1s with positive selection signal are more likely
# to be in vs out of genes
fisher.test(L1SingletonCoeffs$Dist2Gene == 0, L1SingletonCoeffs$blnSelect)

cat("done!\n")

##########################################
#                                        #
#   Regress intersection with genes      #
#                                        #
##########################################

cat("Regressing intersection with genes ...")

#######
# Regress against L1 start
#######

# Indicator variable for intersection with genes
L1_1000G$blnOLGene  <- overlapsAny(L1_1000G_GR_hg19, GeneGR)
L1_1000G$blnOLProm  <- overlapsAny(L1_1000G_GR_hg19, PromGR)
L1_1000G$blnOLExon  <- overlapsAny(L1_1000G_GR_hg19, ExonGR)
L1_1000G$blnOLCpG   <- overlapsAny(L1_1000G_GR_hg19, CpGGR)
L1_1000G$Dist2CpG   <- Dist2Closest(L1_1000G_GR_hg19, CpGGR)
L1_1000G$L1StartNum <- as.numeric(as.character(L1_1000G$L1Start))
L1_1000G$L1EndNum   <- as.numeric(as.character(L1_1000G$L1End))
L1_1000G$blnFull    <- L1_1000G$L1StartNum <= 1 & L1_1000G$L1EndNum >= 6000
hist(L1_1000G$L1StartNum, breaks = 0:6100)
hist(L1_1000G$L1StartNum)
hist(L1_1000G$L1EndNum, breaks = 0:6100)
hist(L1_1000G$InsLength, breaks = 0:6100)
hist(L1_1000G$InsLength)
hist(width(L1RefGR), breaks = 0:6500)
max(width(L1_1000G_GR_hg19))
sum(L1_1000G$L1StartNum == 2, na.rm = T)
sum(L1_1000G$blnOLProm)
mean(L1_1000G$blnOLGene)
table(L1_1000G$blnOLCpG, L1_1000G$blnFull)
t.test(L1_1000G$Dist2CpG ~ L1_1000G$blnFull)
wilcox.test(L1_1000G$Dist2CpG ~ L1_1000G$blnFull)


# Create transparent point color
PCol  <- rgb(0,0,0, alpha = 0.1)
PCol1 <- rgb(1,0,0, alpha = 0.1)
PCol2 <- rgb(0,0,1, alpha = 0.1)

# Regress indicator for gene overlap against L1 start abd frequency
L1_1000G_subset <- L1_1000G[L1_1000G$Frequency >= 0.01, ]
L1_1000G_subset <- L1_1000G
InsLorder       <- order(L1_1000G_subset$L1StartNum)
LogRegGeneOLVsL1Start <- glm(blnOLGene ~ L1StartNum + Frequency + blnFull, 
                             family = "binomial", data = L1_1000G_subset)
summary(LogRegGeneOLVsL1Start)
LogRegGeneOLVsL1Start_noFreq <- glm(blnOLGene ~ L1StartNum, 
                                    family = "binomial", 
                                    data = L1_1000G_subset)
summary(LogRegGeneOLVsL1Start_noFreq)

# Regress indicator for promoter overlap against L1 start abd frequency
LogRegPromOLVsL1Start <- glm(blnOLProm ~ L1StartNum + Frequency, 
                             family = "binomial", data = L1_1000G)
summary(LogRegPromOLVsL1Start)
sum(L1_1000G$blnOLProm)

# Regress combined indicators 
LogRegCombinedOLVsL1Start <- glm((blnOLGene | blnOLProm) ~ L1StartNum*blnOLGene + 
                                   Frequency*blnOLGene, 
                             family = "binomial", data = L1_1000G)
summary(LogRegCombinedOLVsL1Start)

# Proportion of L1 overlapping with genes per L1 start
L1_1000G$L1StartBins <- cut(L1_1000G$L1StartNum, breaks = seq(1, 6021, 20))
PropOLperL1Start <- aggregate(blnOLGene ~ L1StartBins, data = L1_1000G, 
                              FUN = mean)
NperL1Start      <- aggregate(blnOLGene ~ L1StartBins, data = L1_1000G, 
                              FUN = length)
StartperL1Start  <- aggregate(L1StartNum ~ L1StartBins, data = L1_1000G, 
                              FUN = mean)
colnames(NperL1Start)[2] <- "N"
OLperL1Start <- merge(PropOLperL1Start, NperL1Start)
OLperL1Start <- merge(OLperL1Start, StartperL1Start)
OLperL1Start$StDev <- sqrt(OLperL1Start$blnOLGene * (1 - OLperL1Start$blnOLGene) / 
                             OLperL1Start$N)
plot(OLperL1Start$L1StartNum, OLperL1Start$blnOLGene, xlim = c(0, 200))

# Plot proportion gene overlap against L1 start
L1OLVsL1StartSmoothed <- supsmu(L1_1000G_subset$L1StartNum,  
                                1*L1_1000G_subset$blnOLGene)
plot(L1OLVsL1StartSmoothed$x, L1OLVsL1StartSmoothed$y, type = "l", 
     xlab = "L1 5' start",
     ylab = "Proportion of L1s in genes")
lines(L1_1000G_subset$L1StartNum[InsLorder], 
      LogRegGeneOLVsL1Start_noFreq$fitted.values[InsLorder], lty = 2)
CreateDisplayPdf('D:/L1polymORF/Figures/PropInGeneVsInsLength.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 5, width = 5)
jpeg(file = "D:/L1polymORF/Figures/PropInGeneVsInsLength.jpeg")
plot(L1OLVsL1StartSmoothed$x, L1OLVsL1StartSmoothed$y, type = "l", 
     xlab = "L1 5' start",
     ylab = "Proportion of L1s in genes")
lines(L1_1000G_subset$L1StartNum[InsLorder], 
      LogRegGeneOLVsL1Start_noFreq$fitted.values[InsLorder], lty = 2)
dev.off()

plot(L1OLVsL1StartSmoothed$x, L1OLVsL1StartSmoothed$y, type = "l", xlab = "L1 start",
     ylab = "Proportion of L1s in genes", xlim = c(0, 30),
     ylim = c(0.38, 0.55))
# points(L1_1000G$L1StartNum[L1_1000G$blnOLGene], 
#        rep(0.42, sum(L1_1000G$blnOLGene)), col = PCol1, pch = 15, cex = 1.3)
# points(L1_1000G$L1StartNum[!L1_1000G$blnOLGene], 
#        rep(0.41, sum(!L1_1000G$blnOLGene)), col = PCol2, pch = 15, cex = 1.3)
lines(L1_1000G_subset$L1StartNum[InsLorder], 
      LogRegGeneOLVsL1Start_noFreq$fitted.values[InsLorder], col = "red")

# Plot proportion gene overlap against L1 frequency
L1OLVsFreqSmoothed <- supsmu(L1_1000G$Frequency,  1*L1_1000G$blnOLGene)
plot(L1OLVsFreqSmoothed$x, L1OLVsFreqSmoothed$y, type = "l")
points(L1_1000G$Frequency[L1_1000G$blnOLGene], rep(0.4, sum(L1_1000G$blnOLGene)))
plot(L1OLVsFreqSmoothed$x, L1OLVsFreqSmoothed$y, type = "l",
     xlim = c(0, 0.01))
points(L1_1000G$Frequency[L1_1000G$blnOLGene], rep(0.4, sum(L1_1000G$blnOLGene)), col = PCol, 
       pch = 16)

# Frequency and start range with high proportion in genes
blnLowFreq  <- L1_1000G$Frequency >= 0.001 & L1_1000G$Frequency <= 0.003
blnLowStart <- L1_1000G$L1StartNum >= 10 & L1_1000G$L1StartNum <= 15
sum(blnLowFreq & blnLowStart, na.rm = T)
mean(L1_1000G$blnOLGene[blnLowFreq & blnLowStart], na.rm = T)
mean(L1_1000G$blnOLGene[blnLowFreq], na.rm = T)
mean(L1_1000G$blnOLGene[blnLowStart], na.rm = T)

# Regress indicator for exon overlap against L1 start abd frequency
LogRegExonOLVsL1Start <- glm(blnOLExon ~ L1StartNum  + Frequency, 
                             family = "binomial", data = L1_1000G)
summary(LogRegExonOLVsL1Start)

# Regress indicator for exon overlap against L1 start and frequency
# among all L1s in genes
LogRegExonOLVsL1StartFreq <- glm(blnOLExon ~ L1StartNum  + Frequency, 
                             family = "binomial", data = L1_1000G,
                             subset = blnOLGene)
summary(LogRegExonOLVsL1StartFreq)
LogRegExonOLVsL1Start <- glm(blnOLExon ~ L1StartNum, 
                             family = "binomial", data = L1_1000G,
                             subset = blnOLGene)
summary(LogRegExonOLVsL1Start)
LogRegExonOLVsL1Start <- glm(blnOLExon ~ L1StartNum, 
                             family = "binomial", data = L1_1000G,
                             subset = blnOLGene)
summary(LogRegExonOLVsL1Start)

cat("done!\n")

##########################################
#                                        #
#   Regress strand alignedness           #
#                                        #
##########################################

cat("Regressing alignedness ...")

# Overlap map between genes and L1
L1_Gene_OL <- findOverlaps(L1_1000G_GR_hg19, GeneGR)

# Create data.frame that contains L1 start and strandedness for overlapping
# pairs
L1GeneOLInfo <- data.frame(L1Start = L1_1000G$L1StartNum[L1_Gene_OL@from],
                           L1End = L1_1000G$L1EndNum[L1_Gene_OL@from],
                           L1Strand = L1_1000G$L1Strand[L1_Gene_OL@from],
                           GeneStrand = as.vector(strand(GeneGR))[L1_Gene_OL@to])
L1GeneOLInfo$blnSameStrand <- L1GeneOLInfo$L1Strand == L1GeneOLInfo$GeneStrand

# Regress indicator for promoter overlap against L1 start abd frequency
LogRegSameStrandVsL1Start <- glm(blnSameStrand ~ L1Start, 
                             family = "binomial", data = L1GeneOLInfo)
summary(LogRegSameStrandVsL1Start)

# Plot proportion of L1 with same steand as gene against L1 start
L1SameStrandVsL1StartSmoothed <- supsmu(L1GeneOLInfo$L1Start,  1*L1GeneOLInfo$blnSameStrand)
plot(L1SameStrandVsL1StartSmoothed$x, L1SameStrandVsL1StartSmoothed$y, type = "l", xlab = "L1 start",
     ylab = "Proportion of L1s with same strand as gene")
points(L1GeneOLInfo$L1Start[which(L1GeneOLInfo$blnSameStrand)], 
       rep(0.36, sum(L1GeneOLInfo$blnSameStrand, na.rm = T)), col = PCol, pch = 16)
plot(L1SameStrandVsL1StartSmoothed$x, L1SameStrandVsL1StartSmoothed$y, type = "l", xlab = "L1 start",
     ylab = "Proportion of L1s with same strand as gene", xlim = c(0, 50))
points(L1GeneOLInfo$L1Start[which(L1GeneOLInfo$blnSameStrand)], 
       rep(0.36, sum(L1GeneOLInfo$blnSameStrand, na.rm = T)), col = PCol, pch = 16)


##########################################
#                                        #
#   Export IDs of intersecting genes     #
#                                        #
##########################################

# Define a cut-off value for L1 allele frequency
FreqCutOff <- 0.005
blnLowFreq <- L1_1000G$Frequency <= FreqCutOff

# Get genes that overlap with low and high frequency L1a
GeneGR_All      <- subsetByOverlaps(GeneGR, L1_1000G_GR_hg19)
GeneGR_LowFreq  <- subsetByOverlaps(GeneGR, L1_1000G_GR_hg19[blnLowFreq])
GeneGR_HiFreq   <- subsetByOverlaps(GeneGR, L1_1000G_GR_hg19[!blnLowFreq])
GeneGR_Ref      <- subsetByOverlaps(GeneGR, L1RefGR)
PromGR_All      <- subsetByOverlaps(PromGR, L1_1000G_GR_hg19)

PromIDs <- select(Homo.sapiens, 
       keys = as.character(PromGR_All@elementMetadata@listData$tx_id), 
       columns ="ENTREZID", keytype = "TXID")
Prom_ENTREZID <- PromIDs$ENTREZID[!is.na(PromIDs$ENTREZID)]
Prom_ENTREZID <- unique(as.character(Prom_ENTREZID))

# Check whether gene size determines overlap
GeneOLCount     <- countOverlaps(GeneGR, L1_1000G_GR_hg19)
GeneOLCount_ref <- countOverlaps(GeneGR, L1RefGR)
GeneOLNegSelectCount <- countOverlaps(GeneGR, 
                          L1SingletonCoeffs_GR[L1SingletonCoeffs$blnNotSelect])
GeneOLPosSelectCount <- countOverlaps(GeneGR, 
                                   L1SingletonCoeffs_GR[L1SingletonCoeffs$blnSelect])
table(GeneOLPosSelectCount)
boxplot(log(width(GeneGR)) ~ GeneOLCount)
plot(width(GeneGR), GeneOLCount, xlab = "Gene length", xlim = c(0, 2*10^6),
     col = PCol, pch = 16)

sort(width(GeneGR), decreasing = T)[1:5]
max(1*GeneOLCount)

# Write out gene ids for different genelists
writeLines(GeneGR_All@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_All")
writeLines(GeneGR_Ref@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_Ref")
writeLines(GeneGR_LowFreq@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_LowFreq")
writeLines(GeneGR_HiFreq@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_HiFreq")
writeLines(Prom_ENTREZID, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_Promoter")

###################################################
#                                                 #
#   Sample genes at random according to width     #
#                                                 #
###################################################

# Sampled gene indices
idxSampledGenes <- sample(length(GeneGR), 
                         sum(L1_1000G$blnOLGene, na.rm = T), prob = width(GeneGR))
SampledGeneGR <- GeneGR[idxSampledGenes]
writeLines(SampledGeneGR@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/SampledGeneIDs")

cat("done!\n")

###################################################
#                                                 #
#   Analyze enrichment of annotation terms        #
#                                                 #
###################################################

cat("Analyzing enrichment of annotation terms ...")
# Map genomic ranges of gene expression data to gene ranges
GeneGexpOL       <- findOverlaps(GExpGR, GeneGR)
GExpByGeneTissue <- aggregate.data.frame(
  GExpByTissue[GeneGexpOL@from, GexpTissue$name], 
  by = list(GeneGexpOL@to), FUN = mean)

# Create table with gene names
keytypes(org.Hs.eg.db)
cols <- c("UNIPROT")
select(org.Hs.eg.db,
       keys = GeneGR@elementMetadata@listData$gene_id[1:5],
       columns = "IPI", keytype = "ENTREZID")
GeneLookup <- select(org.Hs.eg.db,
                    keys = GeneGR@elementMetadata@listData$gene_id,
                    columns = "UNIPROT", keytype = "ENTREZID",
                    multiVals = "first")
GeneIDmatch <- match(GeneGR@elementMetadata@listData$gene_id, GeneLookup$ENTREZID)
GeneLookup1 <- GeneLookup[GeneIDmatch,]
  
# Get indicator for UniProt IDs
idxUniProt <- which(!is.na(GeneLookup1$UNIPROT))

# Create a data.frame with UniProt ID, gene length and indicator for keyword 
# Membrane
GeneTable <- data.frame(UniProtID = GeneLookup1$UNIPROT[idxUniProt],
                        GeneLength = log10(width(GeneGR)[idxUniProt]),
                        GeneLength_untrans = width(GeneGR)[idxUniProt],
                        L1Count = GeneOLCount[idxUniProt],
                        L1Count_ref = GeneOLCount_ref[idxUniProt],
                        L1NegSelectCount = GeneOLNegSelectCount[idxUniProt],
                        L1PosSelectCount = GeneOLPosSelectCount[idxUniProt],
                        NrMappedGene = NrMappedGene[idxUniProt],
                        NrMappedRange = NrMappedRange[idxUniProt],
                        idxGeneGR = idxUniProt)
GeneTable$blnMHC <- overlapsAny(GeneGR[idxUniProt], MHC_GR)
# Create a uniprot object for humans
IDMatch   <- match(GeneTable$UniProtID, UniProtData$Entry)
GeneTable <- GeneTable[!is.na(IDMatch), ]
IDMatch   <- IDMatch[!is.na(IDMatch)]

# Create indicator variable for each enriched term
ImmGlob <- as.character(InterProData$ENTRY_AC[InterProData$ENTRY_NAME == 
                                                "Immunoglobulin I-set"])
GeneTable$Membrane <- 1:nrow(GeneTable) %in% 
                          grep("Membrane", UniProtData$Keywords[IDMatch])
GeneTable$HostReceptor <- 1:nrow(GeneTable) %in% 
  grep("Host cell receptor for virus entry", UniProtData$Keywords[IDMatch])
GeneTable$Glycoprotein <- 1:nrow(GeneTable) %in% 
  grep("Glycoprotein", UniProtData$Keywords[IDMatch])
GeneTable$CellJunction <- 1:nrow(GeneTable) %in% 
  grep("Cell junction", UniProtData$Keywords[IDMatch])
GeneTable$Kinase <- 1:nrow(GeneTable) %in% 
  grep("Kinase", UniProtData$Keywords[IDMatch])
GeneTable$ImmGlob <- 1:nrow(GeneTable) %in% 
  grep(ImmGlob, UniProtData$Cross.reference..InterPro.[IDMatch])
GeneTable$Fibronect <- 1:nrow(GeneTable) %in% 
  grep("IPR003961", UniProtData$Cross.reference..InterPro.[IDMatch])
GeneTable$Galactose <- 1:nrow(GeneTable) %in% 
  grep("IPR008979", UniProtData$Cross.reference..InterPro.[IDMatch])
GeneTable$Pleckstrin <- 1:nrow(GeneTable) %in% 
  grep("IPR011993", UniProtData$Cross.reference..InterPro.[IDMatch])
GeneTable$Cytokin <- 1:nrow(GeneTable) %in% 
  grep("IPR026791", UniProtData$Cross.reference..InterPro.[IDMatch])
GeneTable$PDZ <- 1:nrow(GeneTable) %in% 
  grep("IPR001478", UniProtData$Cross.reference..InterPro.[IDMatch])
GeneTable$Concanavalin <- 1:nrow(GeneTable) %in% 
  grep("IPR013320", UniProtData$Cross.reference..InterPro.[IDMatch])
GeneTable$Ligand <- 1:nrow(GeneTable) %in% 
  grep("IPR001828", UniProtData$Cross.reference..InterPro.[IDMatch])
GeneTable$MAM <- 1:nrow(GeneTable) %in% 
  grep("IPR000998", UniProtData$Cross.reference..InterPro.[IDMatch])
GeneTable$PTP <- 1:nrow(GeneTable) %in% 
  grep("IPR000242", UniProtData$Cross.reference..InterPro.[IDMatch])
GeneTable$PDEase <- 1:nrow(GeneTable) %in% 
  grep("IPR002073", UniProtData$Cross.reference..InterPro.[IDMatch])

# Perform poisson regression to determine whether enrichment terms affect L1 
# insertion
GLM_OLCount <- glm(L1Count ~ GeneLength + Membrane + Glycoprotein + CellJunction +
                   Kinase + ImmGlob + Fibronect + Galactose + Pleckstrin + Cytokin +
                   PDZ + Concanavalin + Ligand + MAM + PTP + PDEase + HostReceptor +
                     NrMappedGene + NrMappedRange, 
                   data = GeneTable,
                   family = "poisson")
Sum_GLM_OLCount <- summary(GLM_OLCount)
Padj <- p.adjust(Sum_GLM_OLCount$coefficients[,'Pr(>|z|)'])
cbind(Sum_GLM_OLCount$coefficients[Padj < 0.05, 'Estimate'], Padj[Padj < 0.05])
write.csv(cbind(Sum_GLM_OLCount$coefficients, Padj),
          "D:/L1polymORF/Data/AnnoTermPValues.csv")

# Perform poisson regression to determine whether enrichment terms affect L1 
# reference insertion
GLM_OLCount_ref <- glm(L1Count_ref ~ GeneLength + Membrane + Glycoprotein + CellJunction +
                     Kinase + ImmGlob + Fibronect + Galactose + Pleckstrin + Cytokin +
                     PDZ + Concanavalin + Ligand + MAM + PTP + PDEase + HostReceptor, 
                   data = GeneTable,
                   family = "poisson")
Sum_GLM_OLCount_ref <- summary(GLM_OLCount_ref)
Padj_ref <- p.adjust(Sum_GLM_OLCount_ref$coefficients[,'Pr(>|z|)'])
cbind(Sum_GLM_OLCount_ref$coefficients[Padj_ref < 0.05, 'Estimate'], 
      Padj_ref[Padj_ref < 0.05])

# Test whether insertions into glycoproteins are enriched among L1 insertions
# with negative selection when compared to all genic L1
L1_1000G$blnNegSel <- overlapsAny(L1_1000G_GR_hg19, 
                                  L1SingletonCoeffs_GR[L1SingletonCoeffs$blnNotSelect])
table(L1_1000G$blnNegSel, L1_1000G$blnOLProm)

# Test whether glycoproteins are significantly enriched among negatively 
# selected L1
fisher.test(GeneTable$Glycoprotein, GeneTable$L1NegSelectCount)

# Check info on gene with positive selection signal
GeneTable[GeneTable$L1PosSelectCount == 1, ]

# Create 
GLM_OLCount_L <- glm(L1Count ~ GeneLength,  data = GeneTable,
                   family = "poisson")

# Plot smoothed proportion overlap against gene length
OLVsGeneLSmoothed <- supsmu(width(GeneGR),  GeneOLCount)
OLVsGeneLSmoothed_log10 <- supsmu(log10(width(GeneGR)),  GeneOLCount)
plot(log10(width(GeneGR)),  GeneOLCount, col = PCol, pch = 16,
     xaxt = "n",
     xlab = "Gene length [bp]", 
     ylab = "Number of L1 insertions per gene")

axis(1, at = 2:7, paste("10^", c(2:7), sep = ""))
GeneLOrder <- order(GeneTable$GeneLength)
lines(GeneTable$GeneLength[GeneLOrder], 
      GLM_OLCount_L$fitted.values[GeneLOrder], col = "red")
CreateDisplayPdf('D:/L1polymORF/Figures/NrInsVsGeneLength.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 5, width = 5)
jpeg(file = "D:/L1polymORF/Figures/NrInsVsGeneLength.jpeg")
plot(log10(width(GeneGR)),  GeneOLCount, col = PCol, pch = 16,
     xaxt = "n",
     xlab = "Gene length [bp]", 
     ylab = "Number of L1 insertions per gene")

axis(1, at = 2:7, paste("10^", c(2:7), sep = ""))
GeneLOrder <- order(GeneTable$GeneLength)
lines(GeneTable$GeneLength[GeneLOrder], 
      GLM_OLCount_L$fitted.values[GeneLOrder], col = "red")
dev.off()

plot(width(GeneGR),  GeneOLCount, col = PCol, pch = 16,
     xlab = "Gene length [bp]", 
     ylab = "Number of L1 insertions per gene", xlim = c(0, 3*10^6))

axis(1, at = 2:7, paste("10^", c(2:7), sep = ""))
lines(GeneTable$GeneLength_untrans[GeneLOrder], 
      GLM_OLCount_L$fitted.values[GeneLOrder], col = "red")
CreateDisplayPdf('D:/L1polymORF/Figures/NrInsVsGeneLength_untrans.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 5, width = 5)

# Append gene expression data
idxGexpMatch <- match(GeneTable$idxGeneGR, GExpByGeneTissue$Group.1)
GeneTable    <- cbind(GeneTable, GExpByGeneTissue[idxGexpMatch, ])

# Perform poisson regression to test for 
GLMExpr <- paste('glm(L1Count ~ GeneLength + Membrane + Glycoprotein + CellJunction +',
                 'Kinase + ImmGlob + Fibronect + Galactose + Pleckstrin + Cytokin +',
                 'PDZ + Concanavalin + Ligand + MAM + PTP + PDEase + HostReceptor +', 
                 paste(GexpTissue$name, collapse = " + "),
                 ', data =  GeneTable, family = "poisson")')
GLM_OLCount <- eval(parse(text = GLMExpr))
Sum_GLM_OLCount <- summary(GLM_OLCount)
Padj <- p.adjust(Sum_GLM_OLCount$coefficients[,'Pr(>|z|)'])
Sum_GLM_OLCount$coefficients[Padj < 0.05, 'Estimate']

# Compare the fit of log-transformed vs untransformed gene length
max(GeneTable$GeneLength)
GLM_OLCount1 <- glm(L1Count ~ GeneLength, data = GeneTable,
                    family = "poisson", subset = GeneLength <= 6)
GLM_OLCount2 <- glm(L1Count ~ GeneLength_untrans, data = GeneTable,
                    family = "poisson", subset = GeneLength <= 6)
summary(GLM_OLCount1)
summary(GLM_OLCount2)
cat("done!\n")

###################################################
#                                                 #
#   Analyze multiple L1 per gene      #
#                                                 #
###################################################

###################################################
#                                                 #
#   Regress L1 overlap against gene expression    #
#                                                 #
###################################################

GExpByTissue$L1Count <- countOverlaps(GExpGR, L1_1000G_GR_hg19)
GExpByTissue$GeneWidth <- width(GExpGR)

GLMExpr <- paste('glm(L1Count ~ GeneWidth +', paste(GexpTissue$name, collapse = " + "),
                ', data =  GExpByTissue, family = "poisson")')

GLM_L1_Exp <- eval(parse(text = GLMExpr))
summary(GLM_L1_Exp)
GLM_L1_Exp$effects
Sum_GLM_L1_Exp <- summary(GLM_L1_Exp)
Padj <- p.adjust(Sum_GLM_L1_Exp$coefficients[,'Pr(>|z|)'])
Padj[Padj < 0.05]
Sum_GLM_L1_Exp$coefficients[Padj < 0.05, 'Estimate']

###################################################
#                                                 #
#         Analyze glycoproteins                   #
#                                                 #
###################################################

# Create a boolean vector for glycoproteins with L1 insertion
blnL1Glyco <- GeneTable$Glycoprotein & GeneTable$L1Count
sum(blnL1Glyco)

# Get gene range indices of glycoproteins with L1 insertion
idxGlyco   <- GeneTable$idxGeneGR[GeneTable$Glycoprotein]
idxGlycoL1 <- GeneTable$idxGeneGR[blnL1Glyco]
length(idxGlycoL1) / length(idxGlyco)
mean(GeneTable$L1Count > 0)
mean(GeneTable$GeneLength)
mean(GeneTable$GeneLength[GeneTable$L1Count > 0])
mean(GeneTable$GeneLength[GeneTable$Glycoprotein])
mean(GeneTable$GeneLength[GeneTable$Glycoprotein & (GeneTable$L1Count > 0)])

# Get an indicator for overlap of glycoproteins
L1_1000G$blnOLGlyco <- overlapsAny(L1_1000G_GR_hg19, GeneGR[idxGlyco])
L1SingletonCoeffs$blnOLGlyco <- overlapsAny(L1SingletonCoeffs_GR, 
                                            GeneGR[idxGlyco])

# Test whether glycoproteins are overrepresented among negatively selected genes
fisher.test(L1_1000G$blnOLGlyco[L1_1000G$blnOLGene], 
            L1_1000G$blnNegSel[L1_1000G$blnOLGene])

# Test whether selected
wilcox.test(Freq ~ blnNotSelect, data = L1SingletonCoeffs)
t.test(Freq ~ blnNotSelect, data = L1SingletonCoeffs)

LML1NegSelVsFreq <- glm(blnNotSelect ~ Freq, data = L1SingletonCoeffs, 
                        weights = 1/se.coef., family = "quasibinomial",
                        subset = blnOLGlyco)
summary(LML1NegSelVsFreq)
LML1NegSelVsFreq <- glm(blnNotSelect ~ Freq, data = L1SingletonCoeffs, 
                        weights = 1/se.coef., family = "quasibinomial")
summary(LML1NegSelVsFreq)

# Regress glycoprotein indicator against frequency and L1 start
LogRegGlycoOLVsL1Start <- glm(blnOLGlyco ~ L1StartNum + Frequency + blnFull, 
                             family = "binomial", data = L1_1000G)
summary(LogRegGlycoOLVsL1Start)
LogRegGlycoOLVsL1StartOnly <- glm(blnOLGlyco ~ L1StartNum, 
                              family = "binomial", data = L1_1000G)
summary(LogRegGlycoOLVsL1StartOnly)
LogRegGlycoOLVsFreqOnly <- glm(blnOLGlyco ~ Frequency, 
                                  family = "binomial", data = L1_1000G)
summary(LogRegGlycoOLVsFreqOnly)


L1GlycoOLVsL1StartSmoothed <- supsmu(L1_1000G$L1StartNum,  1*L1_1000G$blnOLGlyco)
plot(L1GlycoOLVsL1StartSmoothed$x, L1GlycoOLVsL1StartSmoothed$y, type = "l", 
     xlab = "L1 start", ylab = "Proportion of L1s in glycoproteins")
L1_1000G_noGlyco <- L1_1000G[!L1_1000G$blnOLGlyco, ]
L1noGlycoOLVsL1StartSmoothed <- supsmu(L1_1000G_noGlyco$L1StartNum,  1*L1_1000G_noGlyco$blnOLGene)
plot(L1noGlycoOLVsL1StartSmoothed$x, L1noGlycoOLVsL1StartSmoothed$y, type = "l", 
     xlab = "L1 start", ylab = "Proportion of L1s in genes")

# Check whether frequency or start differs between L1 in glycoprotein genes and other genes
L1_inGenes <- L1_1000G[L1_1000G$blnOLGene, ]
t.test(Frequency ~ blnOLGlyco, data = L1_inGenes)
t.test(L1StartNum ~ blnOLGlyco, data = L1_inGenes)
wilcox.test(L1StartNum ~ blnOLGlyco, data = L1_inGenes)
hist(L1_1000G$L1StartNum)

# Regress glycoprotein indicator against frequency and L1 start
LogRegGlycoOLVsL1Start <- glm(blnOLGlyco ~ L1StartNum + Frequency + blnFull, 
                              family = "binomial", data = L1_inGenes)
summary(LogRegGlycoOLVsL1Start)
LogRegGlycoOLVsL1StartOnly <- glm(blnOLGlyco ~ L1StartNum, 
                                  family = "binomial", data = L1_inGenes)
summary(LogRegGlycoOLVsL1StartOnly)
LogRegGlycoOLVsFreqOnly <- glm(blnOLGlyco ~ Frequency, 
                               family = "binomial", data = L1_inGenes)
summary(LogRegGlycoOLVsFreqOnly)

# Get entrez IDs of glycoproteins 
GeneGR[idxGlyco]
keytypes(org.Hs.eg.db)
GycoProtL1 <- select(org.Hs.eg.db,
       keys = GeneGR@elementMetadata@listData$gene_id[idxGlycoL1],
       columns = c("SYMBOL", "GENENAME"), keytype = "ENTREZID")
write.csv(GycoProtL1, file = "D:/L1polymORF/Data/GlycoProtsWithL1.csv")


###################################################
#                                                 #
#     Predict L1 insertion based on DNAse         #
#                                                 #
###################################################

# Function to create genomic ranges for a particular chromosome number
CreateGR <- function(i, RangeWidth){
  CL <- ChromLengthsHg19[i]
  GRStarts <- seq(1, CL, RangeWidth)
  GREnds  <- c(GRStarts[-1] - 1, CL)
  GRanges(seqnames = names(CL), IRanges(start = GRStarts, end = GREnds))
}

# Create genomic ranges for DNAse and other summaries
SummaryGR <- CreateGR(1, RangeWidth)
for (i in 2:length(ChromLengthsHg19)){
  NewGR <- CreateGR(i, RangeWidth)
  SummaryGR <- c(SummaryGR, NewGR)
}
warnings()
# Summarize score by stem-cells and non-stem cells (move to separate script)
# ScoreByStem <- t(sapply(1:nrow(DNAseData), function(i){
#   IDs <- strsplit(as.character(DNAseData$sourceIds[i]), ",")[[1]]
#   Scores <- as.numeric(strsplit(as.character(DNAseData$sourceScores[i]), ",")[[1]])
#   c(StemScores = sum(Scores[IDs %in% StemIDs]),
#     NotStemScores = sum(Scores[IDs %in% NotStemIDs]))
# }))
dim(ScoreByStem)


# Create genomic ranges for the DNAse data and then aggregate by SummaryGR
DNAseGR <- makeGRangesFromDataFrame(DNAseData)
DNAseSummaryOL <- findOverlaps(DNAseGR, SummaryGR)
hist(width(DNAseGR), breaks = seq(0, 9000, 20), xlim = c(0, 1000))

# Create a factor that indicates whether a DNAse peak occurs in stem cells only
# in non-stem cells only, in both, or neither
blnStem    <- ScoreByStem[,"StemScores"] > 0
blnNotStem <- ScoreByStem[,"NotStemScores"] > 0
DNAseCell  <- rep("neither", nrow(ScoreByStem))
DNAseCell[blnStem & blnNotStem]    <- "both"
DNAseCell[blnStem & (!blnNotStem)] <- "stem"
DNAseCell[(!blnStem) & blnNotStem] <- "diff"
table(blnStem, blnNotStem)

# Get overlap to L1
blnOLDNAseL1 <- overlapsAny(DNAseGR, L1_1000G_GR_hg19)
L1DNAseCombo <- table(blnOLDNAseL1, DNAseCell)
barplot(L1DNAseCombo[2, ] / colSums(L1DNAseCombo))
table(DNAseCell)
chisq.test(blnOLDNAseL1, DNAseCell)

# Aggregate DNAse scores by overlap with L1
DNAseAgg <- aggregate.data.frame(ScoreByStem, by = list(blnOLDNAseL1),
                                 FUN = mean)
DNAseAgg$StemScores / rowSums(DNAseAgg[,2:3])

# Calculate mean DNAse per summary genomic range
DNAseBySummaryGR <- aggregate.data.frame(ScoreByStem[DNAseSummaryOL@from, ], 
                                         by = list(DNAseSummaryOL@to), FUN = mean)
DNAseSummary <- as.data.frame(matrix(0, nrow = length(SummaryGR), ncol = 2))
DNAseSummary[DNAseBySummaryGR$Group.1, ] <- DNAseBySummaryGR[,2:3]
colnames(DNAseSummary) <- colnames(DNAseBySummaryGR)[2:3]

# Add other variables to summary
DNAseSummary$L1Count   <- countOverlaps(SummaryGR, L1_1000G_GR_hg19)
SummaryViews <- BSgenomeViews(BSgenome.Hsapiens.UCSC.hg19, SummaryGR)
GCFreq <- letterFrequency(SummaryViews, letters = "GC")
DNAseSummary$GC <- GCFreq
DNAseSummary$L1Count_fragm   <- countOverlaps(SummaryGR, 
                                   L1_1000G_GR_hg19[which(L1_1000G$InsLength < 5900)])
DNAseSummary$L1Count_full   <- countOverlaps(SummaryGR, 
                                             L1_1000G_GR_hg19[which(L1_1000G$InsLength >= 6000)])
DNAseSummary$ScoreSum <- DNAseSummary$StemScores + DNAseSummary$NotStemScores
DNAseSummary$pStem    <- DNAseSummary$StemScores / (DNAseSummary$ScoreSum + 1)
cor(DNAseSummary$ScoreSum, DNAseSummary$pStem)
cor(DNAseSummary$StemScores, DNAseSummary$NotStemScores)

hist(width(L1_1000G_GR_hg19@ranges))
hist(DNAseSummary$L1Count_full)
hist(DNAseSummary$L1Count_fragm)
# GLM_L1_DNAse_sum <- glm(L1Count ~ ScoreSum + pStem, data = DNAseSummary,
#                          family = "poisson")
# summary(GLM_L1_DNAse_sum)

# Regress LINE-1 count against DNAse scores
GLM_L1_DNAse_both <- glm(L1Count ~ NotStemScores + StemScores + GC, 
                         data = DNAseSummary, family = poisson)
SUM_GLM_L1_DNAse_both <- summary(GLM_L1_DNAse_both)
SUM_GLM_L1_DNAse_both$coefficients[,'Pr(>|z|)']
exp(SUM_GLM_L1_DNAse_both$coefficients[,'Estimate'])
plot(DNAseSummary$NotStemScores, DNAseSummary$L1Count, col = PCol)
plot(DNAseSummary$NotStemScores, DNAseSummary$L1Count, col = PCol)

# Regress LINE-1 count against DNAse scores
ChrV <- as.vector(seqnames(SummaryGR))
ChrV[ChrV == "chrX"] <- "chr23"
ChrV[ChrV == "chrY"] <- "chr24"
DNAseSummary$ChrNum <- as.numeric(substr(ChrV, 4, nchar(ChrV)))
GLM_L1_DNAse_gee <- gee(L1Count ~ NotStemScores + StemScores + GC, 
                         data = DNAseSummary, family = "poisson",
                         corstr =  "exchangeable", Mv = 1, 
                        id = ChrNum)
SUM_GLM_L1_DNAse_gee <- summary(GLM_L1_DNAse_gee)
coef(SUM_GLM_L1_DNAse_gee)[,'Estimate']
se <- SUM_GLM_L1_DNAse_gee$coefficients[, "Robust S.E."]
Ps <- 1 - pnorm(abs(coef(SUM_GLM_L1_DNAse_gee)[,'Estimate']) / se)
format(Ps, digits = 22)
cbind(coef(SUM_GLM_L1_DNAse_gee)[,'Estimate'] - se * qnorm(0.9995),
      coef(SUM_GLM_L1_DNAse_gee)[,'Estimate'] + se * qnorm(0.9995))

GLM_L1_DNAse_stem <- glm(L1Count ~ StemScores, data = DNAseSummary,
                    family = "poisson")
summary(GLM_L1_DNAse_stem)
GLM_L1_DNAse_NotStem <- glm(L1Count ~ NotStemScores, data = DNAseSummary,
                    family = "poisson")
GLM_L1_DNAse_GC <- glm(L1Count ~ GC, data = DNAseSummary,
                            family = "poisson")
summary(GLM_L1_DNAse_GC)
GLM_L1_DNAse_fragm <- glm(L1Count_fragm ~ StemScores + NotStemScores, data = DNAseSummary,
                    family = "poisson")
summary(GLM_L1_DNAse_fragm)
GLM_L1_DNAse_full <- glm(L1Count_full ~ StemScores + NotStemScores, data = DNAseSummary,
                          family = "poisson")
summary(GLM_L1_DNAse_full)
max(DNAseSummary$L1Count_full)

plot(DNAseSummary$StemScores, DNAseSummary$NotStemScores, col = PCol)
plot(DNAseSummary$StemScores, DNAseSummary$L1Count, col = PCol)
plot(DNAseSummary$NotStemScores, DNAseSummary$L1Count, col = PCol)

cor(DNAseSummary$StemScores, DNAseSummary$NotStemScores)
cor(DNAseSummary$StemScores, DNAseSummary$GC)
cor(DNAseSummary$NotStemScores, DNAseSummary$GC)


# Determine proportion of L1 overlapping with heterochromatin
# blnOLHetero <- overlapsAny(L1_1000G_GR_hg19, HeteroGR)
# mean(blnOLHetero)
# sum(width(HeteroGR)/10^6) / sum(ChromLengthsHg19/10^6)
# 
# fisher.test(blnOLHetero, L1_1000G$InsLength >= 6000)

##############################################
#                                            #
#   Analyze difference between               #
#   expected and observed heterozygosity     #
#                                            #
##############################################

# Get expected and observed number of homozygotes
L1_1000G$HomoExp <- length(SampleColumns) * L1_1000G$Frequency^2
L1_1000G$HomoObs <- rowSums(L1_1000G[,SampleColumns] == 2)
L1_1000G$HomoDiff <- (L1_1000G$HomoObs - L1_1000G$HomoExp) 
mean(L1_1000G$HomoDiff, na.rm = T)
hist(L1_1000G$HomoDiff, breaks = -500:200)
hist(L1_1000G$HomoExp, breaks = 0:2000, xlim = c(0, 10)) 

# Check whether differemce between expected and observed homozygosity differes
# between L1 in genes and outside genes
boxplot(HomoDiff ~ blnOLGene, data = L1_1000G)
t.test(HomoDiff ~ blnOLGene, data = L1_1000G)
wilcox.test(HomoDiff ~ blnOLGene, data = L1_1000G)
with(L1_1000G, mean(HomoDiff[blnOLGene], na.rm = T))
with(L1_1000G, mean(HomoDiff[!blnOLGene], na.rm = T))

##########################################
#                                        #
#   Test gene-singleton association      #
#                                        #
##########################################

# Create boolean variable for intersection with genes
L1SingletonCoeffs$blnOLGene <- L1SingletonCoeffs$Dist2Gene == 0
fisher.test(L1SingletonCoeffs$blnOLGene, L1SingletonCoeffs$blnNotSelect)
wilcox.test(CoefSt ~ blnOLGene, data = L1SingletonCoeffs)
t.test(CoefSt ~ blnOLGene, data = L1SingletonCoeffs)
aggregate(CoefSt ~ blnOLGene, data = L1SingletonCoeffs, FUN = mean)

wilcox.test(coef ~ blnOLGene, data = L1SingletonCoeffs)
t.test(coef ~ blnOLGene, data = L1SingletonCoeffs)
aggregate(coef ~ blnOLGene, data = L1SingletonCoeffs, FUN = mean)

sum(L1_1000G$blnOLGene) / 3060

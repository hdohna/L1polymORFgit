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
GeneGRPath      <- "D:/L1polymORF/Data/Gene_GR.RData"

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
load(ParalogPath)
load(GeneGRPath)

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
ExonGR <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
PromGR <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream = 10000)
CDSGR  <- cds(TxDb.Hsapiens.UCSC.hg19.knownGene)
IntronGRList   <- intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                      use.names = T)
FiveUTRGRList  <- fiveUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                       use.names = T)
ThreeUTRGRList <- threeUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                        use.names = T)
# Among overlapping genomic ranges, retain the longest
GeneGR <- GeneGR_NCBIRefSeq

# Read UniProt data
UniProtData  <- read.delim(UniProtPath)
UniProtData  <- UniProtData[,1:5]
InterProData <- read.delim(InterProPath)

# Read in table with L1HS from the reference genome
L1RefTab <- read.csv(L1RefPath)
L1RefGR <- makeGRangesFromDataFrame(L1RefTab, seqnames.field = "genoName",
                                    start.field = "genoStart",
                                    end.field = "genoEnd")

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

LML1NegSelVsL1StartFreq <- glm(blnNegSelect ~ Freq + L1Start, 
                           data = L1SinglCoeff_nonzeroSE, 
                           weights = 1/se.coef., family = "quasibinomial")
summary(LML1NegSelVsL1StartFreq)
LML1NegSelVsL1Start <- glm(blnNegSelect ~ L1Start + Freq, 
                           data = L1SinglCoeff_nonzeroSE, 
                      weights = 1/se.coef., family = "quasibinomial")
summary(LML1NegSelVsL1Start)

LML1NegSelVsFreq <- glm(blnNegSelect ~ Freq, data = L1SinglCoeff_nonzeroSE, 
                   weights = 1/se.coef., family = "quasibinomial")
summary(LML1NegSelVsFreq)

# Determine whether L1s with negative selection signal are more likely
# to be in vs out of genes
table(L1SingletonCoeffs$Dist2Gene == 0, L1SingletonCoeffs$blnNegSelect)
fisher.test(L1SingletonCoeffs$Dist2Gene == 0, L1SingletonCoeffs$blnNegSelect)
chisq.test(L1SingletonCoeffs$Dist2Gene == 0, L1SingletonCoeffs$blnNegSelect)

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

# Indicator variable for intersection with various GRanges
L1_1000G$blnOLGene  <- overlapsAny(L1_1000G_GR_hg19, GeneGR, ignore.strand = T)
L1_1000G$blnOLGeneSameStrand <- overlapsAny(L1_1000G_GR_hg19, GeneGR)
L1_1000G$blnOLProm  <- overlapsAny(L1_1000G_GR_hg19, PromGR, ignore.strand = T)
L1_1000G$blnOLExon  <- overlapsAny(L1_1000G_GR_hg19, ExonGR, ignore.strand = T)
L1_1000G$L1StartNum <- as.numeric(as.character(L1_1000G$L1Start))
L1_1000G$L1EndNum   <- as.numeric(as.character(L1_1000G$L1End))
L1_1000G$blnFull    <- L1_1000G$L1StartNum <= 1 & L1_1000G$L1EndNum >= 6000

# Overlap map between genes and L1
L1_Gene_OL <- findOverlaps(L1_1000G_GR_hg19, GeneGR, ignore.strand = T)

# Add Gene info
L1_1000G$GeneWidth <- NA
L1_1000G$GeneWidth[L1_Gene_OL@from] <- width(GeneGR)[L1_Gene_OL@to]

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

# Regress indicator for promoter overlap against L1 start abd frequency
LogRegPromOLVsL1Start <- glm(blnOLProm ~ L1StartNum + Frequency, 
                             family = "binomial", data = L1_1000G)
summary(LogRegPromOLVsL1Start)
sum(L1_1000G$blnOLProm)

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
LogRegGeneOLVsL1Start_noFreq <- glm(blnOLGene ~ L1StartNum, 
                                    family = "binomial", data = L1_1000G)
L1OLVsL1StartSmoothed <- supsmu(L1_1000G_subset$L1StartNum,  
                                1*L1_1000G_subset$blnOLGene)
plot(L1OLVsL1StartSmoothed$x, L1OLVsL1StartSmoothed$y, type = "l", 
     xlab = "L1 5' start",
     ylab = "Proportion of L1s in genes")
lines(L1_1000G_subset$L1StartNum[InsLorder], 
      LogRegGeneOLVsL1Start_noFreq$fitted.values[InsLorder], lty = 2)
CreateDisplayPdf('D:/L1polymORF/Figures/PropInGeneVsInsLength_Gencode.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 5, width = 5)
jpeg(file = "D:/L1polymORF/Figures/PropInGeneVsInsLength_Gencode.jpeg")
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

# Frequency and start range with high proportion in genes
blnLowFreq  <- L1_1000G$Frequency >= 0.001 & L1_1000G$Frequency <= 0.003
blnLowStart <- L1_1000G$L1StartNum >= 10 & L1_1000G$L1StartNum <= 15
sum(blnLowFreq & blnLowStart, na.rm = T)
mean(L1_1000G$blnOLGene[blnLowFreq & blnLowStart], na.rm = T)
mean(L1_1000G$blnOLGene[blnLowFreq], na.rm = T)
mean(L1_1000G$blnOLGene[blnLowStart], na.rm = T)


cat("done!\n")

##########################################
#                                        #
#   Regress strand alignedness           #
#                                        #
##########################################

cat("Regressing alignedness ...")

# Create data.frame that contains L1 start and strandedness for overlapping
# pairs
blnNoDupl <- !duplicated(L1_Gene_OL@from)
L1GeneOLInfo <- data.frame(L1Start = L1_1000G$L1StartNum[L1_Gene_OL@from[blnNoDupl]],
                           L1End = L1_1000G$L1EndNum[L1_Gene_OL@from[blnNoDupl]],
                           L1Strand = L1_1000G$L1Strand[L1_Gene_OL@from[blnNoDupl]],
                           Freq = L1_1000G$Frequency[L1_Gene_OL@from[blnNoDupl]],
                           GeneStrand = as.vector(strand(GeneGR))[L1_Gene_OL@to[blnNoDupl]],
                           GeneWidth = width(GeneGR)[L1_Gene_OL@to[blnNoDupl]])
L1GeneOLInfo$blnSameStrand <- L1GeneOLInfo$L1Strand == L1GeneOLInfo$GeneStrand
sum(is.na(L1GeneOLInfo$blnSameStrand))


# Test whether alignedness differs significantly from 0.5
blnNotNA <- !is.na(L1GeneOLInfo$blnSameStrand)
NrSame   <- sum(L1GeneOLInfo$blnSameStrand[blnNotNA])
pbinom(q = NrSame, size = sum(blnNotNA), prob = 0.5)
NrSame / sum(blnNotNA)

##########################################
#                                        #
#   Export IDs of intersecting genes     #
#                                        #
##########################################

# Define a cut-off value for L1 allele frequency
FreqCutOff <- 0.005
blnLowFreq <- L1_1000G$Frequency <= FreqCutOff

# Get genes that overlap with low and high frequency L1a
GeneGR_All        <- subsetByOverlaps(GeneGR, L1_1000G_GR_hg19, ignore.strand = T)
GeneGR_SameStrand <- subsetByOverlaps(GeneGR, L1_1000G_GR_hg19)
GeneGR_LowFreq    <- subsetByOverlaps(GeneGR, L1_1000G_GR_hg19[blnLowFreq], ignore.strand = T)
GeneGR_HiFreq     <- subsetByOverlaps(GeneGR, L1_1000G_GR_hg19[!blnLowFreq], ignore.strand = T)
GeneGR_Ref        <- subsetByOverlaps(GeneGR, L1RefGR, ignore.strand = T)
PromGR_All        <- subsetByOverlaps(PromGR, L1_1000G_GR_hg19, ignore.strand = T)

PromIDs <- select(Homo.sapiens, 
       keys = as.character(PromGR_All@elementMetadata@listData$tx_id), 
       columns ="ENTREZID", keytype = "TXID")
Prom_ENTREZID <- PromIDs$ENTREZID[!is.na(PromIDs$ENTREZID)]
Prom_ENTREZID <- unique(as.character(Prom_ENTREZID))

# Check whether gene size determines overlap
GeneOLCount     <- countOverlaps(GeneGR, L1_1000G_GR_hg19, ignore.strand = T)
GeneOLCount_ref <- countOverlaps(GeneGR, L1RefGR, ignore.strand = T)
GeneOLNegSelectCount <- countOverlaps(GeneGR, 
                          L1SingletonCoeffs_GR[L1SingletonCoeffs$blnNegSelect],
                          ignore.strand = T)
GeneOLPosSelectCount <- countOverlaps(GeneGR, 
                                   L1SingletonCoeffs_GR[L1SingletonCoeffs$blnSelect],
                                   ignore.strand = T)
table(GeneOLPosSelectCount)
boxplot(log(width(GeneGR)) ~ GeneOLCount)
plot(width(GeneGR), GeneOLCount, xlab = "Gene length", xlim = c(0, 2*10^6),
     col = PCol, pch = 16)

sort(width(GeneGR), decreasing = T)[1:5]
max(1*GeneOLCount)

# Write out gene ids for different genelists
writeLines(as.character(GeneGR_All@elementMetadata@listData$name2),
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_All_Gencode")
writeLines(as.character(GeneGR_SameStrand@elementMetadata@listData$name2),
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_SameStrand_Gencode")


###################################################
#                                                 #
#   Analyze enrichment of annotation terms        #
#                                                 #
###################################################

cat("Analyzing enrichment of annotation terms ...")

# Create a gene id lookup table
GeneLookup <- select(org.Hs.eg.db,
                     keys = as.character(GeneGR@elementMetadata@listData$name2),
                     columns = "UNIPROT", keytype = "SYMBOL",
                     multiVals = "first")
GeneIDmatch <- match(GeneGR@elementMetadata@listData$name2, GeneLookup$SYMBOL)
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
                        idxGeneGR = idxUniProt)
GeneTable$blnMHC <- overlapsAny(GeneGR[idxUniProt], MHC_GR, ignore.strand = T)


# Create a uniprot object for humans
IDMatch   <- match(GeneTable$UniProtID, UniProtData$Entry)
GeneTable <- GeneTable[!is.na(IDMatch), ]
IDMatch   <- IDMatch[!is.na(IDMatch)]

# Create indicator variables for enriched Uniprot keywords
GeneTable$Membrane <- 1:nrow(GeneTable) %in% 
                          grep("Membrane", UniProtData$Keywords[IDMatch])
GeneTable$HostReceptor <- 1:nrow(GeneTable) %in% 
  grep("Host cell receptor for virus entry", UniProtData$Keywords[IDMatch])
GeneTable$AltSplice <- 1:nrow(GeneTable) %in% 
  grep("Alternative splicing", UniProtData$Keywords[IDMatch])
GeneTable$Polymorph <- 1:nrow(GeneTable) %in% 
  grep("Polymorphism", UniProtData$Keywords[IDMatch])
GeneTable$Membrane <- 1:nrow(GeneTable) %in% 
  grep("Membrane", UniProtData$Keywords[IDMatch])
GeneTable$Glycoprotein <- 1:nrow(GeneTable) %in% 
  grep("Glycoprotein", UniProtData$Keywords[IDMatch])
GeneTable$CellJunction <- 1:nrow(GeneTable) %in% 
  grep("Cell junction", UniProtData$Keywords[IDMatch])
GeneTable$Kinase <- 1:nrow(GeneTable) %in% 
  grep("Kinase", UniProtData$Keywords[IDMatch])
GeneTable$Synapse <- 1:nrow(GeneTable) %in% 
  grep("Synapse", UniProtData$Keywords[IDMatch])
GeneTable$Transmembrane <- 1:nrow(GeneTable) %in% 
  grep("Transmembrane", UniProtData$Keywords[IDMatch])
GeneTable$CellAdh <- 1:nrow(GeneTable) %in% 
  grep("Cell adhesion", UniProtData$Keywords[IDMatch])


# Create indicator variables for enriched Interpro terms
ImmGlob <- as.character(InterProData$ENTRY_AC[InterProData$ENTRY_NAME == 
                                                "Immunoglobulin I-set"])
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
colnames(GeneTable)
sum(GeneTable$blnMHC)
GLM_OLCount <- glm(L1Count ~ GeneLength + AltSplice + Polymorph + Membrane + Glycoprotein + CellJunction +
                   Kinase + ImmGlob + Fibronect + Galactose + Pleckstrin + Cytokin +
                   PDZ + Concanavalin + Ligand + MAM + PTP + PDEase + HostReceptor +
                   blnMHC, 
                   data = GeneTable,
                   family = "poisson")
Sum_GLM_OLCount <- summary(GLM_OLCount)
Padj <- p.adjust(Sum_GLM_OLCount$coefficients[,'Pr(>|z|)'])
cbind(Sum_GLM_OLCount$coefficients[Padj < 0.05, 'Estimate'], Padj[Padj < 0.05])
write.csv(cbind(Sum_GLM_OLCount$coefficients, Padj),
          "D:/L1polymORF/Data/AnnoTermPValues_NCBIRef.csv")

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
                                  L1SingletonCoeffs_GR[L1SingletonCoeffs$blnNegSelect],
                                  ignore.strand = T)
table(L1_1000G$blnNegSel, L1_1000G$blnOLProm)

# Test whether glycoproteins are significantly enriched among negatively 
# selected L1
fisher.test(GeneTable$Glycoprotein, GeneTable$L1NegSelectCount)
fisher.test(GeneTable$Glycoprotein & GeneTable$Membrane, GeneTable$L1NegSelectCount)

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
CreateDisplayPdf('D:/L1polymORF/Figures/NrInsVsGeneLength_NCBIRef.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 5, width = 5)
jpeg(file = "D:/L1polymORF/Figures/NrInsVsGeneLength_NCBIRef.jpeg")
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
CreateDisplayPdf('D:/L1polymORF/Figures/NrInsVsGeneLength_untrans_NCBIRef.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 5, width = 5)

###################################################
#                                                 #
#         Analyze glycoproteins                   #
#                                                 #
###################################################

cat("Analyze glycoproteins ...")

# Create a boolean vector for glycoproteins with L1 insertion
blnL1Glyco <- GeneTable$Glycoprotein & GeneTable$L1Count
sum(blnL1Glyco)

# Get gene range indices of glycoproteins with L1 insertion
idxGlyco   <- GeneTable$idxGeneGR[GeneTable$Glycoprotein]
idxGlycoMembrane  <- GeneTable$idxGeneGR[GeneTable$Glycoprotein & GeneTable$Membrane]
idxGlycoL1 <- GeneTable$idxGeneGR[blnL1Glyco]
length(idxGlycoL1) / length(idxGlyco)
mean(GeneTable$L1Count > 0)
mean(GeneTable$GeneLength)
mean(GeneTable$GeneLength[GeneTable$L1Count > 0])
mean(GeneTable$GeneLength[GeneTable$Glycoprotein])
mean(GeneTable$GeneLength[GeneTable$Glycoprotein & (GeneTable$L1Count > 0)])

# Get an indicator for overlap of glycoproteins
L1_1000G$blnOLGlyco <- overlapsAny(L1_1000G_GR_hg19, GeneGR[idxGlyco], 
                                   ignore.strand = T)
L1_1000G$blnOLGlycoMem <- overlapsAny(L1_1000G_GR_hg19, GeneGR[idxGlycoMembrane], 
                                   ignore.strand = T)
L1SingletonCoeffs$blnOLGlyco <- overlapsAny(L1SingletonCoeffs_GR, 
                                          GeneGR[idxGlyco], ignore.strand = T)
L1SingletonCoeffs$blnOLGlycoMem <- overlapsAny(L1SingletonCoeffs_GR, 
                                            GeneGR[idxGlycoMembrane], ignore.strand = T)

# Test whether glycoproteins are overrepresented among negatively selected genes
fisher.test(L1_1000G$blnOLGlyco[L1_1000G$blnOLGene], 
            L1_1000G$blnNegSel[L1_1000G$blnOLGene])
fisher.test(L1_1000G$blnOLGlycoMem[L1_1000G$blnOLGene], 
            L1_1000G$blnNegSel[L1_1000G$blnOLGene])

# Test whether glycoproteins are overrepresented among genes that have the same
# strand as L1
blnStrandInfo <- !is.na(L1_1000G$L1Strand)
blnMaxWidth   <- L1_1000G$GeneWidth <= 5*10^5
fisher.test(L1_1000G$blnOLGlyco[L1_1000G$blnOLGene & blnStrandInfo], 
            L1_1000G$blnOLGeneSameStrand[L1_1000G$blnOLGene & blnStrandInfo])
fisher.test(L1_1000G$blnOLGlycoMem[L1_1000G$blnOLGene & blnStrandInfo], 
            L1_1000G$blnOLGeneSameStrand[L1_1000G$blnOLGene & blnStrandInfo])
table(L1_1000G$blnOLGlyco[L1_1000G$blnOLGene & blnStrandInfo], 
            L1_1000G$blnOLGeneSameStrand[L1_1000G$blnOLGene & blnStrandInfo])
fisher.test(L1_1000G$blnOLGlycoMem[which(L1_1000G$blnOLGene & blnStrandInfo & blnMaxWidth)], 
            L1_1000G$blnOLGeneSameStrand[which(L1_1000G$blnOLGene & blnStrandInfo & blnMaxWidth)])

LML1NegSelVsFreq_all <- glm(blnNegSelect ~ Freq, data = L1SingletonCoeffs, 
                             weights = 1/se.coef., family = "quasibinomial")
summary(LML1NegSelVsFreq_all)
LML1NegSelVsFreq_glyc <- glm(blnNegSelect ~ Freq, data = L1SingletonCoeffs, 
                        weights = 1/se.coef., family = "quasibinomial",
                        subset = blnOLGlycoMem)
summary(LML1NegSelVsFreq_glyc)
LML1NegSelVsFreq <- glm(blnNegSelect ~ Freq, data = L1SingletonCoeffs, 
                        weights = 1/se.coef., family = "quasibinomial")
summary(LML1NegSelVsFreq)

# Regress glycoprotein indicator against frequency and L1 start
LogRegGlycoOLVsL1Start <- glm(blnOLGlyco ~ L1StartNum + Frequency + blnFull, 
                             family = "binomial", data = L1_1000G,
                             subset = blnOLGene)
summary(LogRegGlycoOLVsL1Start)
LogRegGlycoOLVsL1StartOnly <- glm(blnOLGlyco ~ L1StartNum, 
                              family = "binomial", data = L1_1000G, 
                              subset = blnOLGene)
summary(LogRegGlycoOLVsL1StartOnly)
LogRegGlycoOLVsFreqOnly <- glm(blnOLGlyco ~ Frequency, 
                                  family = "binomial", data = L1_1000G,
                               subset = blnOLGene)
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
wilcox.test(Frequency ~ blnOLGlyco, data = L1_inGenes)
wilcox.test(L1StartNum ~ blnOLGlyco, data = L1_inGenes)

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
cat("done!\n")
cat("Analysis completed!\n")


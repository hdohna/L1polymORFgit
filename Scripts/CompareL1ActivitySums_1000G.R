# The script below reads GRanges with L1s from the 1000 Genome data
# Phase 3, available at ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
# The script takes matches L1 insertions to 
# The Granges were created by the script 'Create_1000G_L1GRanges.R'

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(DescTools)
library(ade4)
library(vegan)
library(ggplot2)
library(grid)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Specify file paths
G1000SamplePath <- 'D:/L1polymORF/Data/1000GenomeSampleInfo.txt'
L1GRPath        <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
L1RefRangePath  <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'

# Specify parameters
NrInfoCols <- 9

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

# Load previously generated objects
load(L1RefRangePath)
load(L1GRPath)

##########
# Process L1 catalog
##########

cat("Matching 1000G L1 to L1 catalog\n")

# Read in table with Info about 1000 genome samples 
SampleInfo_1000Genome <- read.table(G1000SamplePath, as.is = T, header = T)
table(SampleInfo_1000Genome$super_pop)

# Read in table with known L1 
L1Catalogue <- read.csv("D:/L1polymORF/Data/L1CatalogExtended.csv", as.is = T)
L1Catalogue$Allele[is.na(L1Catalogue$Allele)] <- 1

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
blnInRef1         <- (L1Catalogue$end_HG38 - L1Catalogue$start_HG38) > 6000 
L1CatalogL1Mapped <- L1Catalogue[blnL1Mapped & blnAllele1,]

# Create genomic ranges for catalog L1
L1CatalogGR <- GRanges(seqnames = L1CatalogL1Mapped$Chromosome,
                       ranges = IRanges(start = pmin(L1CatalogL1Mapped$start_HG38,
                                                     L1CatalogL1Mapped$end_HG38),
                                        end = pmax(L1CatalogL1Mapped$start_HG38,
                                                   L1CatalogL1Mapped$end_HG38)),
                       strand = L1CatalogL1Mapped$strand_L1toRef)

DistCat2_1000G <- Dist2Closest(L1CatalogGR, L1_1000G_GRList_hg38$LiftedRanges)

# Get indices of 1000 Genome and catalog elements that match
idx1000G <- nearest(L1CatalogGR, L1_1000G_GRList_hg38$LiftedRanges)
L1CatalogMatch1000G <- L1CatalogL1Mapped[DistCat2_1000G < 100, ]
idx1000GMatchCat    <- idx1000G[DistCat2_1000G < 100]
L1CatalogGRMatched  <- L1CatalogGR[DistCat2_1000G < 100]
sum(DistCat2_1000G < 100)
sum(DistCat2_1000G < 10)

# add a column with dummy activity sums
L1CatalogMatch1000G$ActivityDummy <- 1

# Get rows of L1_1000G table that can be matched to catalog
L1_1000G_match <- L1_1000G[idx1000GMatchCat, ]

# Get expected and observed activity sums
ExpectedAct <- 2 * L1CatalogMatch1000G$ActivityNum %*% L1_1000G_match$Frequency
ObservedAct <- sapply(SampleColumns, function(x){
  L1CatalogMatch1000G$ActivityNum %*% L1_1000G_match[,x]})
hist(ObservedAct)
mean(ObservedAct)
length(ObservedAct)
mean(L1_1000G_reduced$Frequency)

##########################################
#                                        #
#     Get number of insertions around    #
#        matched LINE-1s                 #
#                                        #
##########################################

# Create large genomic ranges around each L1
RangeSeq <- seq(2 * 10^7, 5*10^7, 10^7/2)
L1CorMat <- sapply(RangeSeq, function(Range){
  cat("Calculating correlations for range", Range, "\n")
  L1CatGR_large <- resize(L1CatalogGRMatched, Range, fix = "center")
  OverlapL1Cat  <- findOverlaps(L1CatGR_large, 
                               L1_1000G_GRList_hg38$LiftedRanges[-idx1000GMatchCat])
  
  # Create a vector that counts for each matched L1 the difference between L1 in
  # the vicinity between individuals carrying the L1 and the ones that don't
  idxOverlap <- unique(OverlapL1Cat@from)
  L1Cor <- sapply(idxOverlap, function(i){
    FocalSum   <- L1_1000G_match[i, SampleColumns]
    idxOtherL1 <- OverlapL1Cat@to[OverlapL1Cat@from == i]
    L1other    <- colSums(L1_1000G[idxOtherL1, SampleColumns])
    cor(t(FocalSum), L1other)  
  })        
  
})
ChromLengths
dim(L1CorMat)
plot(RangeSeq, colSums(L1CorMat))
ChromCount <- table(L1CatalogMatch1000G$Chromosome)
SingleL1Chroms <- names(ChromCount)[ChromCount == 1]
idxHot <- which(L1CatalogMatch1000G$ActivityNum > 10 &
                  L1CatalogMatch1000G$Chromosome %in% SingleL1Chroms)
Cols <- rainbow(length(idxHot))
plot(RangeSeq, L1CorMat[idxHot[1], ], type = "l", col = Cols[1], 
     ylim = c(-0.1, 0.5), xlab = "Range width", ylab = "correlation")
for (i in 2:length(idxHot)){
  lines(RangeSeq, L1CorMat[idxHot[i], ], col = Cols[i] )
}
legend("topleft", col = Cols, lty = 1, 
       legend = L1CatalogMatch1000G$Chromosome[idxHot])

# Plot the difference between L1s in samples with and without focal L1 against
# activity of focal L1
# plot(L1CatalogMatch1000G$ActivityNum[idxOverlap], L1Diff)
# cor.test(L1CatalogMatch1000G$ActivityNum[idxOverlap], L1Diff)

#OverlapL1Cat <- findOverlaps()
# OverlapCount <- rep(0, nrow(L1CatalogMatch1000G))
# names(OverlapCount) <- 1:length(OverlapCount)
# OverlapCount_tmp <- table(OverlapL1Cat@from)
# OverlapCount[names(OverlapCount_tmp)] <- OverlapCount_tmp
# plot(L1CatalogMatch1000G$ActivityNum, OverlapCount)
# cor.test(L1CatalogMatch1000G$ActivityNum, OverlapCount)

##########################################
#                                        #
#     Correlation betw. total retortrans #
#     and other L1 insertions            #
#                                        #
##########################################

# The frequency columns
#FreqCols <- (grep("FORMAT", colnames(L1_1000G)) + 1):(ncol(L1_1000G) -3)

# Create a boolean matrix for frequencies above 0
blnFreqAbove0 <- L1_1000G[,SampleColumns] > 0

# Get insertions that occur in only one individual
blnSingl   <- rowSums(blnFreqAbove0) == 1
blnAbove20 <- rowSums(blnFreqAbove0) >= 20
idxAbove20 <- setdiff(which(blnAbove20), idx1000GMatchCat)
idxAll     <- setdiff(1:nrow(L1_1000G), idx1000GMatchCat)

# Get number of single insertion per genome and test for correlation 
# with activity sum
NrSingle <- colSums(blnFreqAbove0[blnSingl, ])
cor.test(ObservedAct, NrSingle, method = "spearman")
plot(ObservedAct, NrSingle)
GLM_NrSingle <- glm(NrSingle ~ ObservedAct, family = poisson)
summary(GLM_NrSingle)

# Get number of insertion per genome and test for correlation 
# with activity sum
NrAbove20 <- colSums(blnFreqAbove0[idxAbove20, ])
cor.test(ObservedAct, NrAbove20, method = "spearman")
plot(ObservedAct, NrAbove20)
GLM_NrAbove20 <- glm(NrAbove20 ~ ObservedAct, family = poisson)
summary(GLM_NrAbove20)

# Get number of insertion per genome and test for correlation 
# with activity sum
NrAll <- colSums(blnFreqAbove0[idxAll, ])
cor.test(ObservedAct, NrAll, method = "spearman")
plot(ObservedAct, NrAll)
GLM_NrAll <- glm(NrAll ~ ObservedAct, family = poisson)
summary(GLM_NrAll)

# Create random activity sums and calculate correlation between 
# activity and insertion frequency
NrSamples <- 500
RandCor <- sapply(1:NrSamples, function(x){
  ActRand <- sample(L1CatalogMatch1000G$ActivityNum, nrow(L1CatalogMatch1000G))
  ObservedActRand <- sapply(SampleColumns, function(x){
    ActRand %*% L1_1000G_match[,x]})
  c(CorAll = cor(ObservedActRand, NrAll, method = "spearman"),
    CorAbove20 = cor(ObservedActRand, NrAbove20, method = "spearman"),
    CorSingle = cor(ObservedActRand, NrSingle, method = "spearman"))
})
sum(RandCor["CorAll", ]     >= cor(ObservedAct, NrAll, method = "spearman"))/NrSamples
sum(RandCor["CorAbove20", ] >= cor(ObservedAct, NrAbove20, method = "spearman"))/NrSamples
sum(RandCor["CorSingle", ]  >= cor(ObservedAct, NrSingle, method = "spearman"))/NrSamples

hist(RandCor["CorAll", ], main = "", xlab = "Correlation",
     breaks = seq(-0.1, 0.5, 0.025))
segments(cor(ObservedAct, NrAll, method = "spearman"), y0 = 0, y1 = 1000, 
         col = "red")
par(mfrow = c(1, 1))
CreateDisplayPdf(Filenamepath = 'D:/L1polymORF/Figures/L1ActSumCorrelations.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

#################################################
#                                               #
#  Correlation betw. L1 activity per chromsome  #
#     and other L1 insertions                   #
#                                               #
#################################################

# Create a data.frame that keeps track of correlation and P-Value
CorPerChr <- data.frame(Chrom = unique(L1CatalogMatch1000G$Chromosome),
                        Cor = NA, P = NA, Slope = NA)
PlotList <- vector(mode = "list", length = nrow(CorPerChr))
for(i in 1:nrow(CorPerChr)){
  ChrName <- as.character(CorPerChr$Chrom[i])
  cat("Analyzing", ChrName, "\n")
  Chr        <- substr(ChrName, 4, nchar(ChrName))
  ChrName    <- paste("chr", Chr, sep = "")
  idxChr     <- which(L1_1000G$CHROM == Chr)
  blnL1Cat   <- L1CatalogMatch1000G$Chromosome == ChrName
  idxMatchL1 <- intersect(idxChr, idx1000GMatchCat)
  idxChr     <- setdiff(idxChr, idx1000GMatchCat)
  L1Act <- sapply(SampleColumns, function(x){
    L1CatalogMatch1000G$ActivityNum[blnL1Cat] %*% L1_1000G_match[blnL1Cat,x]
  })
  OtherL1 <- colSums(L1_1000G[idxChr,SampleColumns])
  PlotData <- data.frame(L1Act = L1Act, OtherL1 = OtherL1)
  LM <- lm(OtherL1 ~ L1Act, data = PlotData)
  
  PlotList [[i]] <- ggplot(PlotData, aes(x=L1Act, y=OtherL1)) + 
    geom_point(alpha = 0.05, size = 3) +
    geom_abline(slope = LM$coefficients[2], intercept = LM$coefficients[1]) +
    scale_x_continuous(name = "") + 
    scale_y_continuous(name = "") +
    ggtitle(ChrName)
  
  CT <- cor.test(L1Act, OtherL1)
  CorPerChr$Cor[i] <- CT$estimate
  CorPerChr$P[i] <- CT$p.value
  CorPerChr$Slope[i] <- LM$coefficients[2]
}

multiplot(plotlist = PlotList[!is.na(CorPerChr$Cor)], cols = 4)
hist(CorPerChr$P)

# Create a matrix of random correlations per chromosome
NrSamples <- 500
table(L1CatalogMatch1000G$Chromosome)
Chroms2Analyze <- CorPerChr$Chrom[!is.na(CorPerChr$Cor)]
ChrName <- "chr2"
RandCorsPerChrom <- sapply(as.character(Chroms2Analyze), function(ChrName){
  cat("Getting random correlations for", ChrName, "\n")
  Chr <- substr(ChrName, 4, nchar(ChrName))
  idxChr     <- which(L1_1000G$CHROM == Chr)
  blnL1Cat   <- L1CatalogMatch1000G$Chromosome == ChrName
  RandCor <- sapply(1:NrSamples, function(x){
    idxMatchL1 <- sample(idxChr, sum(blnL1Cat))
    idxOther   <- setdiff(idxChr, idxMatchL1)
    OtherL1    <- colSums(L1_1000G[idxOther, SampleColumns])
    Activity   <- L1CatalogMatch1000G$ActivityNum[blnL1Cat]
    ObservedActRand <- sapply(SampleColumns, function(x){
      Activity %*% L1_1000G[idxMatchL1, x]})
    cor(ObservedActRand, OtherL1)
  })
})
sum(RandCor_chr2 >= CorPerChr$Cor[2]) / NrSamples

# Match chromosomes
ChrMatch   <- match(Chroms2Analyze, CorPerChr$Chrom)
blnGreater <- t(RandCorsPerChrom) >= CorPerChr$Cor[ChrMatch]
rowSums(blnGreater) / NrSamples
hist(rowSums(blnGreater) / NrSamples, breaks = seq(0, 1, 0.05))
sum(rowMeans(RandCorsPerChrom, na.rm = T) >= mean(CorPerChr$Cor[ChrMatch])) / 
  NrSamples

# Matrix of correlation between L1 frequencies
CM <- cor(t(L1_1000G_match[, SampleColumns]))
hist(CM)
sum(CM == 1, na.rm = T)
diag(CM) <- NA
mean(CM, na.rm = T)

# Calculate covariance matrix
CVM <- cov(t(L1_1000G_match[, SampleColumns]))
diag(CVM) <- NA
mean(CVM, na.rm = T)

# Matrix of distances between L1s.
DM <- sapply(1:length(L1CatalogGRMatched), function(i){
  distance(L1CatalogGRMatched, L1CatalogGRMatched[i])
})
DM[is.na(DM)] <- max(DM, na.rm = T)
diag(DM) <- NA
mantel(1 - abs(CM), DM)
plot(DM, CM)

# Matrix of activity sum
AM <- outer(X = L1CatalogMatch1000G$ActivityNum, 
            Y = L1CatalogMatch1000G$ActivityNum,
      FUN = '+')
mantel.partial(1 - CM, AM, DM)
mantel(1 - CM, AM)

#################################################
#                                               #
#       Difference between expected and         #
#             observed L1 counts                #
#                                               #
#################################################

i <- 1
ObsFreqNull <- rep(0, 3)
names(ObsFreqNull) <- 0:2
NTot <- length(SampleColumns)
FreqDiffs <- sapply(1:nrow(L1_1000G_match), function(i){
  ni <- sum(L1_1000G_match[i, SampleColumns])
  FreqV <- unlist(L1_1000G_match[i, SampleColumns])
  FreqsTable <- table(FreqV)
  ObsFreqs   <- ObsFreqNull
  ObsFreqs[names(FreqsTable)] <- FreqsTable
  ExpFreqs <- round(c((NTot - ni)^2, (NTot - ni)*ni, ni^2) / NTot)
  (ObsFreqs -ExpFreqs) %*% c(0:2)
})
hist(FreqDiffs)
L1_1000G_match[,2510:2516]
plot(rowMeans(FreqDiffs))
boxplot(t(FreqDiffs))

AllCounts <- 0:2
Slopes <- sapply(1:ncol(FreqDiffs), function(i){
  LM <- lm(FreqDiffs[,i] ~ AllCounts)
  LM$coefficients[2]
})
plot(L1CatalogMatch1000G$ActivityNum, FreqDiffs)
cor.test(L1CatalogMatch1000G$ActivityNum, FreqDiffs, method = "spearman")

##########################################
#                                        #
#     Compare frequency and 
#     distance to genes among fragments  #
#                                        #
##########################################

# Get ranges of genes
GRgenes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Auxiliary function to get distances to closest gene
DistGeneObj <- distanceToNearest(L1_1000G_GR_hg19, GRgenes_hg19, ignore.strand = T) 
Dist2ClosestGene <- DistGeneObj@elementMetadata@listData$distance

LMF <- lm(L1_1000G_reduced$Frequency ~ Dist2ClosestGene + L1_1000G_reduced$InsLength+
            Dist2ClosestGene : L1_1000G_reduced$InsLength)
summary(LMF)
cor.test(L1_1000G_reduced$Frequency, Dist2ClosestGene, method = "kendall")

##########################################
#                                        #
#     Analyze linkage disequilibrium     #
#     (across all samples)               #
#                                        #
##########################################

# Get linkage disequilibrium, activity sum and indicator of same chromosome 
# for pairs of LINE-1s
cat("Calculating linkage disequilibrium\n")
LDMat <- data.frame(LD = NA, Cor = NA, SameChrom = NA, ActSum = NA, ChiSq = NA)
x <- 1
y <- 2
for(x in 1:(nrow(L1_1000G_match) - 1)){
  blnA  <- L1_1000G_match[x, SampleColumns] > 0
#  Pa    <- L1_1000G_match$Frequency[x]
  Pa    <- mean(blnA)
  for(y in (x + 1):nrow(L1_1000G_match)){
#    Pb           <- L1_1000G_match$Frequency[y]
#    PaPb         <- (2*(Pa * (1 - Pa) + Pa^2)) * (2*(Pb * (1 - Pb) + Pb^2))
#    PaPb         <- Pa*Pb
    blnSameChrom <- L1_1000G_match$CHROM[x] == L1_1000G_match$CHROM[y]
    blnB         <- L1_1000G_match[y, SampleColumns] > 0
    Pb           <- mean(blnB)
    Pab          <- mean(blnA * blnB)
    ChiSq        <- chisq.test(x = as.vector(blnA), y = as.vector(blnB))$statistic
    Cor <- cor(t(L1_1000G_match[x, SampleColumns]), 
               t(L1_1000G_match[y, SampleColumns]))
    if (Cor > 1) browser()
    LDMat   <- rbind(LDMat, c(LD = (Pab - Pa*Pb), 
                              Cor = Cor, 
                              SameChrom = blnSameChrom,
                              ActSum  = L1CatalogMatch1000G$ActivityNum[x] +
                                 L1CatalogMatch1000G$ActivityNum[y],
                              ChiSq = ChiSq))
  }
}
LDMat <- LDMat[!is.na(LDMat$LD),]
max( LDMat$ActSum)
max( LDMat$LD)
LDMat[LDMat$LD == max(LDMat$LD), ]

# Plot histogram of linkage disequilibrium values
hist(LDMat$LD)
boxplot(LD ~ SameChrom, data = LDMat)
boxplot(ActSum ~ SameChrom, data = LDMat)
plot(LD ~ ActSum, data = LDMat)
plot(LDMat$ActSum, LDMat$Cor)
plot(LDMat$LD, LDMat$Cor)
points(LDMat$Cor[LDMat$SameChrom], LDMat$Cor[LDMat$SameChrom])
cor.test(LDMat$LD, LDMat$ActSum, method = "spearman")
cor.test(LDMat$Cor, LDMat$ActSum, method = "spearman")
t.test(LD ~ SameChrom, data = LDMat)
LMLD <- lm(LD ~ SameChrom + ActSum, data = LDMat)
LMCor<- lm(Cor ~ SameChrom + ActSum + SameChrom*ActSum, data = LDMat)
LMCorSubset <- lm(Cor ~ SameChrom + ActSum, data = LDMat,
           subset = Cor < 1)
summary(LMLD)
summary(LMCor)
summary(LMCorSubset)

# Plot activity vs frequency
plot(L1CatalogMatch1000G$ActivityNum, L1_1000G_match$Frequency)
cor(L1CatalogMatch1000G$ActivityNum, L1_1000G_match$Frequency)
cor.test(L1CatalogMatch1000G$ActivityNum, L1_1000G_match$Frequency, 
         method = "spearman")

# Determine difference between expected and observed frequency of
# homozygous L1
L1CatalogMatch1000G$HomoDiff <- sapply(1:nrow(L1_1000G_match), function(x){
  L1_1000G_match$Frequency[x]^2 - mean(L1_1000G_match[x, SampleColumns] == 2)
})
plot(L1CatalogMatch1000G$ActivityNum, L1CatalogMatch1000G$HomoDiff)
cor(L1CatalogMatch1000G$ActivityNum, L1CatalogMatch1000G$HomoDiff)
cor.test(L1CatalogMatch1000G$ActivityNum, L1CatalogMatch1000G$HomoDiff, method = "spearman")

# Compare frequency and activity
plot(L1CatalogMatch1000G$ActivityNum, L1_1000G_match$Frequency)
cor.test(L1CatalogMatch1000G$ActivityNum, L1_1000G_match$Frequency, method = "spearman")

##########################################
#                                        #
#     Sample L1 activity sums            #
#     (across all samples)               #
#                                        #
##########################################


# Sample individual activity sums
cat("Sampling L1 activity sums\n")
NrSamples <- 1000000
SampledActSums <- sapply(1:NrSamples, function(x){
  rbinom(nrow(L1_1000G_match), 2, L1_1000G_match$Frequency) %*% 
    L1CatalogMatch1000G$ActivityNum
})
mean(SampledActSums)

# Get samples of variances
StartVals <- seq(1, NrSamples, length(ObservedAct))
SampledVars <- sapply(StartVals[-length(StartVals)], function(x){
  var(SampledActSums[x:(x+length(ObservedAct))])
})

# Plot variance histogram
par(mar = c(5, 4, 4, 2) + 0.1)
hist(SampledVars, main = "", xlab = "Variance in activity sum")
segments(var(ObservedAct), y0 = 0, y1 = 50, col = "red")
PVar <- mean(SampledVars <= var(ObservedAct))
CreateDisplayPdf('D:/L1polymORF/Figures/L1VarianceOfActivitySums.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')


# Sample quantiles of activity sums
cat("Sampling quantiles of activity sums\n")
SamplQuantList <- SampleQuantiles(SampledActSums, length(SampleColumns), 
                                  NrSamples = 10000, QuantV = seq(0, 1, 0.0001),
                                  LowerQ = 0.005, UpperQ = 0.995)
QSMat <- SamplQuantList$QMat
idxF <- 1:ncol(QSMat)
idxR <- ncol(QSMat):1
mean(SampledActSums)
QQ1  <- qqplot(SampledActSums, ObservedAct, plot.it = F)
plot(QQ1$x, QQ1$y, xlab = "Sampled activity sums", 
     ylab = "Observed activity sums")
MedianMatch1 <- sapply(QQ1$x, function(z) which.min(abs(z - QSMat[1, ])))
MedianMatch2 <- sapply(QQ1$x, function(z) which.min(abs(z + 1 - QSMat[1, ])))
MedianMatch3 <- sapply(QQ1$x, function(z) which.min(abs(z - 1 - QSMat[1, ])))
polygon(QSMat[1, c(idxF, idxR)], c(QSMat[2, ], QSMat[3, idxR]), 
        col = "grey", border = NA)
points(QQ1$x, QQ1$y)
lines(c(0, 1000), c(0, 1000))
# blnOutside <- (QQ1$y < QSMat[2, MedianMatch1] | QQ1$y > QSMat[3, MedianMatch1]) &
#   (QQ1$y < QSMat[2, MedianMatch2] | QQ1$y > QSMat[3, MedianMatch2]) &
#   (QQ1$y < QSMat[2, MedianMatch3] | QQ1$y > QSMat[3, MedianMatch3])
blnOutside <- (QQ1$y < QSMat[2, MedianMatch1] | QQ1$y > QSMat[3, MedianMatch1])
sum(blnOutside)
points(QQ1$x[blnOutside],  QQ1$y[blnOutside], col = "red")
  
CreateDisplayPdf('D:/L1polymORF/Figures/L1ActivitySums.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

# Plot variance histogram and quantiles
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)
hist(SampledVars, xlab = "Variance in activity sum",
     main = "A")
segments(var(ObservedAct), y0 = 0, y1 = 50, col = "red")
PVar <- mean(SampledVars <= var(ObservedAct))

plot(QQ1$x, QQ1$y, xlab = "Sampled activity sums", 
     ylab = "Observed activity sums", main = "B")
MedianMatch1 <- sapply(QQ1$x, function(z) which.min(abs(z - QSMat[1, ])))
MedianMatch2 <- sapply(QQ1$x, function(z) which.min(abs(z + 1 - QSMat[1, ])))
MedianMatch3 <- sapply(QQ1$x, function(z) which.min(abs(z - 1 - QSMat[1, ])))
polygon(QSMat[1, c(idxF, idxR)], c(QSMat[2, ], QSMat[3, idxR]), 
        col = "grey", border = NA)
points(QQ1$x, QQ1$y)
lines(c(0, 1000), c(0, 1000))
blnOutside <- (QQ1$y < QSMat[2, MedianMatch1] | QQ1$y > QSMat[3, MedianMatch1]) &
  (QQ1$y < QSMat[2, MedianMatch2] | QQ1$y > QSMat[3, MedianMatch2]) &
  (QQ1$y < QSMat[2, MedianMatch3] | QQ1$y > QSMat[3, MedianMatch3])
sum(blnOutside)
points(QQ1$x[blnOutside],  QQ1$y[blnOutside], col = "red")

CreateDisplayPdf('D:/L1polymORF/Figures/L1VarianceAndQQActivitySums.pdf', height = 4,
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

# Compare histogram of observed activity sum with histograms of simulated sums
BinWidth <- 25
HistBreaks <- seq(0, 2000, BinWidth)
Obs = ObservedAct
ObsDens <- hist(ObservedAct, HistBreaks, plot = F)$density
SampledDens <- sapply(StartVals[-length(StartVals)], function(x){
  hist(SampledActSums[x:(x+length(ObservedAct))], HistBreaks, plot = F)$density
})
dim(SampledDens)
SampledDensMean <- rowMeans(SampledDens)
SampledDensQMat <- apply(SampledDens, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))
dim(SampledDensQMat)

plot(c(0, HistBreaks[-length(HistBreaks)]), c(0, ObsDens), 
     type = "s", xlim = c(0, 600))
polygon(c(HistBreaks[-length(HistBreaks)], 
          HistBreaks[(length(HistBreaks) - 1):1]), 
        c(SampledDensQMat[1,], SampledDensQMat[3,(length(HistBreaks) - 1):1]), col = "grey",
        border = NA)
lines(c(0, HistBreaks[-length(HistBreaks)]), c(0, SampledDensMean),
      type = "s", col = "red")
lines(c(0, HistBreaks[-length(HistBreaks)]), c(0, ObsDens), 
     type = "s")

plot(HistBreaks[-length(HistBreaks)], (ObsDens - SampledDensMean) / SampledDensMean,
     xlim = c(0, 1000))
lines(c(0, 5000), c(0, 0), col = "red", lty = 2)
mean(ObservedAct)

##########################################
#                                        #
#     Sample L1 activity sums            #
#     (with population structure)         #
#                                        #
##########################################

# Total number of samples
NrSamplesPerPop <- 10000

# Sample individual activity sums
cat("Sampling L1 activity sums with population structure\n")

# Get the number of observed samples per population
NrObsPerPop <- table(SampleInfo_1000Genome$super_pop)

# Get a list of samples per population
SamplePerPopList <-  lapply(names(NrObsPerPop), function(Pop){
  blnPop         <- SampleInfo_1000Genome$super_pop == Pop
  SampleInfo_1000Genome$sample[blnPop]
})
sapply(SamplePerPopList, length)

# Loop over populations and calculate frequency per population
FreqMat <- sapply(SamplePerPopList, function(Samples){
  rowSums(L1_1000G_match[,Samples]) / 2 / length(Samples)
})

# Matrix of sampled L1 activity sums
SampleMat <- sapply(1:NrSamplesPerPop, function(x){
  
  # Initialize sample vector
  SampledActSums_P <- c()
  
  # Loop over populations and sample activity sums
  for (i in 1:length(NrObsPerPop)){
    
    # Sample activity sums for current observation
    SampledActSums <- sapply(1:NrObsPerPop[i], function(x){
      rbinom(nrow(FreqMat), 2, FreqMat[,i]) %*% 
        L1CatalogMatch1000G$ActivityNum
    })
    SampledActSums_P <- c(SampledActSums_P, SampledActSums)
  }
  SampledActSums_P
})
dim(SampleMat)

# Get samples of variances
SampledVars_P <- apply(SampleMat, 2, var)

# Get 99% range of each quantile
QuantV = seq(0, 1, 0.0001) 
LowerQ = 0.005
UpperQ = 0.995
NrSamples = ncol(SampleMat)
SampleMeans <- rep(NA, NrSamples)
SampledQMat <- matrix(nrow = NrSamples, ncol = length(QuantV))
for (i in 1:NrSamples){
    SampleVals <- SampleMat[,i]
    SampleMeans[i]  <- mean(SampleVals)
    SampledQMat[i,] <- quantile(SampleVals, QuantV)
  }
QSMat <- apply(SampledQMat, 2, 
                      FUN = function(x) quantile(x, c(0.5, LowerQ, UpperQ)))
idxF <- 1:ncol(QSMat)
idxR <- ncol(QSMat):1
mean(SampledActSums)
QQ1  <- qqplot(SampleMat[,1], ObservedAct, plot.it = F)
#QQ1  <- qqplot(as.vector(sampleMat), ObservedAct, plot.it = F)

###########
# Plot variance histogram and quantiles
###########

# Set parameters
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)

# Plot variance 
hist(SampledVars_P, main = "A", xlab = "Variance in activity sum",
     xlim = c(10500, 13500))
segments(var(ObservedAct), y0 = 0, y1 = 1000, col = "red")
PVar_P <- mean(SampledVars_P <= var(ObservedAct))

# Plot quantiles
plot(QQ1$x, QQ1$y, xlab = "Sampled activity sums", 
     ylab = "Observed activity sums", main = "B")
MedianMatch1 <- sapply(QQ1$x, function(z) which.min(abs(z - QSMat[1, ])))
MedianMatch2 <- sapply(QQ1$x, function(z) which.min(abs(z + 1 - QSMat[1, ])))
MedianMatch3 <- sapply(QQ1$x, function(z) which.min(abs(z - 1 - QSMat[1, ])))
polygon(QSMat[1, c(idxF, idxR)], c(QSMat[2, ], QSMat[3, idxR]), 
        col = "grey", border = NA)
points(QQ1$x, QQ1$y)
lines(c(0, 1000), c(0, 1000))
blnOutside <- (QQ1$y < QSMat[2, MedianMatch1] | QQ1$y > QSMat[3, MedianMatch1]) &
  (QQ1$y < QSMat[2, MedianMatch2] | QQ1$y > QSMat[3, MedianMatch2]) &
  (QQ1$y < QSMat[2, MedianMatch3] | QQ1$y > QSMat[3, MedianMatch3])
sum(blnOutside)
points(QQ1$x[blnOutside],  QQ1$y[blnOutside], col = "red")

CreateDisplayPdf('D:/L1polymORF/Figures/L1ActivityVarianceAndQQ_PopStr.pdf', 
                 height = 4,
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

# Compare means
mean(colMeans(SampleMat))
mean(ObservedAct)

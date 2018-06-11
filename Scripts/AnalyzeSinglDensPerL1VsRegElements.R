# The script below reads analyzes the singleton density per L1
##########################################
#                                        #
#     Load packages                      #
#                                        #
##########################################

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(survival)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(KernSmooth)
library(glmnet)

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
HiCFolderPath   <- 'D:/L1polymORF/Data/HiCData/'
InputPath       <- 'D:/L1polymORF/Data/SingletonAnalysis.RData'

# Number of info columns in vcf file
NrInfoCols   <- 9

# Minimum number of carriers for a LINE-1 to be analyzed
MinNrCarrier <- 3

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("Loading and processing data ...")

# Load previously generated objects
load(L1GRPath)
load(ChrLPath)
load(InputPath)

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

##########
#  Process ChromHMM
##########

# Read in table with regulatory elements
# RegTable       <- read.table("D:/L1polymORF/Data/ChromHMM", header = T)
# RegTable      <- read.table("D:/L1polymORF/Data/EncodeBroadHMM/wgEncodeBroadHmmH1hescHMM.txt",
#                              col.names = c("bin", "chrom", "ChromStart",
#                                            "ChromEnd", "name", "score", "strand",
#                                            "thickStart", "thickEnd", "itemRgb"))
RegTable      <- read.table("D:/L1polymORF/Data/EncodeBroadHMM/ChromHMMcombined.txt",
                            header = T)
blnAllCellTypes <- RegTable$CellType == 
  "Gm12878,H1hesc,Hepg2,Hmec,Hsmm,Huvec,K562,Nhek,Nhlf"
RegTableAll <- RegTable[blnAllCellTypes,]
#RegTable <- RegTable[blnAllCellTypes,]

idxEnhancer <- grep("Enhancer", RegTable$name)
idxTxn      <- grep("Txn", RegTable$name)
idxProm     <- grep("Promoter", RegTable$name)
idxHetero   <- grep("Heterochrom", RegTable$name)
idxRepr     <- union(grep("Insulator", RegTable$name),
                        grep("Repressed", RegTable$name))

# Counting different regulatory elements
NT <- table(RegTable$name)
sort(NT)

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
RegGR       <- makeGRangesFromDataFrame(RegTable)
RegGRAll    <- makeGRangesFromDataFrame(RegTableAll)
EnhancerGR  <- RegGR[idxEnhancer]
TxnGR       <- RegGR[idxTxn]
PromGR      <- RegGR[idxProm]
ReprGR      <- RegGR[idxRepr]
EnhTxPromGR <- RegGR[unique(c(idxEnhancer, idxTxn, idxProm))]

# check distance of transcription and promoter regions to genes
DistTxn2Gene <- Dist2Closest(TxnGR, 
             genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
sum(DistTxn2Gene == 0)/length(TxnGR)
DistProm2Gene <- Dist2Closest(PromGR, 
                             genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
sum(DistProm2Gene == 0)/length(PromGR)


##########
#  Process transcription factor binding site data
##########

TfbdTable <- read.table("D:/L1polymORF/Data/wgEncodeRegTfbsClusteredV3.txt",
                        col.names = c("bin", "chrom", "start", "end", "name", "score",
                                      "expCount", "expIDs", "expScores"))
TfbGR <- makeGRangesFromDataFrame(TfbdTable)
table(TfbdTable$name)

# Get unique experiment IDs
SplitExpIDList <- lapply(as.character(TfbdTable$expIDs), function(x) as.numeric(strsplit(x, ",")[[1]]))
UniqueExpIDs   <- sort(unique(unlist(SplitExpIDList)))

##########
#  Process vista enhancer data
##########

VistaEnhancerTable <- read.table("D:/L1polymORF/Data/VistaEnhancers.txt", header = T)
VistaEnhancerGR    <- makeGRangesFromDataFrame(VistaEnhancerTable,
                                               start.field = "chromStart",
                                               end.field = "chromEnd")

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
#L1SingletonCoeffs$blnSelect[is.na(L1SingletonCoeffs$blnSelect)] <- FALSE
  
table(L1SingletonCoeffs$blnSelect, L1SingletonCoeffs$blnFull)
fisher.test(L1SingletonCoeffs$blnSelect, L1SingletonCoeffs$blnFull)
nrow(L1SingletonCoeffs)

plot(L1SingletonCoeffs$Freq, L1SingletonCoeffs$coef)
cor.test(L1SingletonCoeffs$Freq[L1SingletonCoeffs$blnSelect],
     L1SingletonCoeffs$coef[L1SingletonCoeffs$blnSelect],
     method = "spearman")

# Indicator for negative selection
L1SingletonCoeffs$blnNotSelect <- L1SingletonCoeffs$blnSig &
  L1SingletonCoeffs$coef > 0

# Indicator fo selection (+1 = positive, -1 = negative, 0 = neutral)
L1SingletonCoeffs$SelectInd <- 0
L1SingletonCoeffs$SelectInd[L1SingletonCoeffs$blnSelect]     <- 1
L1SingletonCoeffs$SelectInd[L1SingletonCoeffs$blnNotSelect]  <- -1

# Bin for insertion length
L1SingletonCoeffs$InsLBins <- cut(L1SingletonCoeffs$InsLength, 
                                  breaks = seq(0, 7000, 1000))
table(L1SingletonCoeffs$InsLBins)

# Caclulate distance to genes
L1SingletonCoeffs$Dist2Gene <- Dist2Closest(L1SingletonCoeffs_GR, 
                                            genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
# Caclulate logarithm of distance to genes
L1SingletonCoeffs$LogDist2Gene <- log(L1SingletonCoeffs$Dist2Gene + 0.1)

# Caclulate distance to closest enhancer
L1SingletonCoeffs$Dist2Enhancer <- Dist2Closest(L1SingletonCoeffs_GR, EnhancerGR)
L1SingletonCoeffs$Dist2EnhancerOrGene <- 
  pmin(L1SingletonCoeffs$Dist2Gene, L1SingletonCoeffs$Dist2Enhancer)

table(L1SingletonCoeffs$Dist2Gene < L1SingletonCoeffs$Dist2Enhancer,
      L1SingletonCoeffs$blnSelect , L1SingletonCoeffs$blnFull)

# Calculate distance to closest enhancer
L1SingletonCoeffs$Dist2Txn <- Dist2Closest(L1SingletonCoeffs_GR, TxnGR)

# Calculate distance to closest promoter
L1SingletonCoeffs$Dist2Prom <- Dist2Closest(L1SingletonCoeffs_GR, PromGR)

# Distance to closest transcription, promoter or enhancer
L1SingletonCoeffs$Dist2EnhTxProm <- Dist2Closest(L1SingletonCoeffs_GR, EnhTxPromGR)

# Distance to closest repressing element
L1SingletonCoeffs$Dist2Repr <- Dist2Closest(L1SingletonCoeffs_GR, ReprGR)
sum(L1SingletonCoeffs$Dist2Repr == 0) / nrow(L1SingletonCoeffs)

(sum((L1SingletonCoeffs$Dist2Repr == 0) & L1SingletonCoeffs$blnFull) / sum(L1SingletonCoeffs$blnFull))^2

# Closest regulatory element
# idxNearestReg <- nearest(L1SingletonCoeffs_GR, ReprGR)
# RegNames <- RegTable$name[-idxHeteroTable]
# L1SingletonCoeffs$ClosestReg <- RegNames[idxNearestReg]

# Calculate distance to closest transcription factor binding site
L1SingletonCoeffs$Dist2Tfb <- Dist2Closest(L1SingletonCoeffs_GR, TfbGR)

# Closest transcription factor binding site
idxNearestTfb  <- nearest(L1SingletonCoeffs_GR, TfbGR)
L1SingletonCoeffs$ClosestTfb <- TfbdTable$name[idxNearestTfb]

L1SingletonCoeffs$ClosestTfb[L1SingletonCoeffs$blnFull & 
                               L1SingletonCoeffs$blnSelect]
L1SingletonCoeffs$ClosestTfb[L1SingletonCoeffs$blnSelect]

TfbRatio <- table(L1SingletonCoeffs$ClosestTfb[L1SingletonCoeffs$blnSelect])/ table(TfbdTable$name)
sort(TfbRatio)

# Caclulate distance to closest vista enhancer
L1SingletonCoeffs$Dist2VistaEnhancer <- Dist2Closest(L1SingletonCoeffs_GR, 
                                                     VistaEnhancerGR)

L1SingletonCoeffs$Dist2VistaOrGene <- pmin(L1SingletonCoeffs$Dist2VistaEnhancer,
                                           L1SingletonCoeffs$Dist2Gene)

##########################################
#                                        #
#        Plot coefficients               #
#                                        #
##########################################

# Plot coefficients vs insertion length
par(mfrow = c(1, 1))
plot(L1SingletonCoeffs$InsLength, L1SingletonCoeffs$coef, 
     xlab = "Insertion length", ylab = "Singleton coefficient")
with(L1SingletonCoeffs, points(InsLength[blnSelect], coef[blnSelect], 
       col = "red"))
MeanCoeffPerL <- aggregate(L1SingletonCoeffs[,c("coef", "InsLength", "blnSelect")],
                           by = list(L1SingletonCoeffs$InsLBins), 
                           FUN = function(x) mean(x, na.rm = T))
plot(MeanCoeffPerL$InsLength, MeanCoeffPerL$coef)
plot(MeanCoeffPerL$InsLength, MeanCoeffPerL$blnSelect)
plot(L1SingletonCoeffs$Freq, L1SingletonCoeffs$coef)
cor.test(L1SingletonCoeffs$Freq, L1SingletonCoeffs$coef)
cor.test(L1SingletonCoeffs$Freq, L1SingletonCoeffs$SelectInd, method = "spearman")
cor.test(L1SingletonCoeffs$Freq, L1SingletonCoeffs$coef)


##########################################
#                                        #
#        Regress coefficients            #
#                                        #
##########################################

#######
# Regress against L1 start
#######

max(L1SingletonCoeffs$L1End, na.rm = T)
# Subset L1 coefficients
sum(is.na(L1SingletonCoeffs$InsLength))
L1SingletonCoeffs_subset <- subset(L1SingletonCoeffs, 
                                   subset = robust.se > 0 & (!is.na(blnSelect))&
                                     L1SingletonCoeffs$L1End >= 5900)
LM_All_Interact_binom <- glm(blnNotSelect ~ L1Start + blnFull  + Freq, 
                             data = L1SingletonCoeffs_subset, 
                             weights = 1/robust.se,
                             family = quasibinomial)
summary(LM_All_Interact_binom)
LM_All_Interact_binom <- glm(1L*blnNotSelect ~ L1Start + blnFull, 
                             data = L1SingletonCoeffs_subset, 
                             weights = 1/robust.se,
                             family = quasibinomial)
summary(LM_All_Interact_binom)


# Smooth proportions of 
blnNotFull <- !L1SingletonCoeffs_subset$blnFull
PropSmoothed_InsL <- supsmu(L1SingletonCoeffs_subset$L1Start[blnNotFull],
                            1*L1SingletonCoeffs_subset$blnSelect[blnNotFull],
                            wt = 1/L1SingletonCoeffs_subset$robust.se[blnNotFull])
PropNSSmoothed_InsL <- supsmu(L1SingletonCoeffs_subset$L1Start[blnNotFull],
                            1*L1SingletonCoeffs_subset$blnNotSelect[blnNotFull],
                            wt = 1/L1SingletonCoeffs_subset$robust.se[blnNotFull])
plot(PropNSSmoothed_InsL$x, PropNSSmoothed_InsL$y, xlab = "L1 insertion length [bp]",
     ylab = "Proportion of L1 with positive selection signal", type = "l",
     col = "red")
lines(PropSmoothed_InsL$x, PropSmoothed_InsL$y)
rect(13, 0, 21, 1, border = "red")
lines(PropSmoothed_InsL$x, PropSmoothed_InsL$y)
InsLorder <- order(L1SingletonCoeffs_subset$L1Start)
lines(L1SingletonCoeffs_subset$L1Start[InsLorder], 
      LM_All_Interact_binom$fitted.values[InsLorder], col = "red")
rect(13, 0, 21, 1, col = "red")
# CreateDisplayPdf('D:/L1polymORF/Figures/PropSelectVsInsLength.pdf', 
#                  PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
#                  height = 5, width = 5)


#######
# Regress against L1 start among rare L1
#######

# Subset L1 coefficients
blnSubset <- L1SingletonCoeffs$robust.se > 0 & 
  (!is.na(L1SingletonCoeffs$blnSelect)) &
  L1SingletonCoeffs$Freq <= 3/2504 &
  L1SingletonCoeffs$L1End >= 5900
L1SingletonCoeffs_subsetLowFreq <- 
  L1SingletonCoeffs[blnSubset, ]
LM_All_Interact_binom_lowFreq <- glm(1L*blnSig ~ L1Start + blnFull, 
                             data = L1SingletonCoeffs_subsetLowFreq, 
                             weights = 1/robust.se,
                             family = quasibinomial)
summary(LM_All_Interact_binom_lowFreq)


#######
# Regress against frequency
#######


LM_SelectFreq_binom <- glm(blnSig ~ Freq + blnSelect:Freq, 
                           data = L1SingletonCoeffs, 
                           weights = 1/robust.se,
                           family = quasibinomial, subset = Freq < 0.02)
LM_NotSelectFreq_binom <- glm(1L*blnNotSelect ~ Freq, 
                              data = L1SingletonCoeffs_subset, 
                              weights = 1/robust.se,
                              family = quasibinomial)
summary(LM_SelectFreq_binom)
summary(LM_NotSelectFreq_binom)
aggregate(Freq ~ SelectInd, data = L1SingletonCoeffs, FUN = mean)

# Plot coefficients against
plot(L1SingletonCoeffs$Freq * 2*length(SampleColumns), 
     L1SingletonCoeffs$coef,
     xlim = c(0, 100))
points(L1SingletonCoeffs$Freq[L1SingletonCoeffs$blnSig] * 2*length(SampleColumns), 
       L1SingletonCoeffs$coef[L1SingletonCoeffs$blnSig], 
       col = "red", pch = 16)



# Plot p-values against frequency
plot(L1SingletonCoeffs$Freq * 2*length(SampleColumns), 
     log(L1SingletonCoeffs$Pr...z..),
     xlim = c(0, 100))
points(L1SingletonCoeffs$Freq[L1SingletonCoeffs$blnSig] * 2*length(SampleColumns), 
       log(L1SingletonCoeffs$Pr...z..[L1SingletonCoeffs$blnSig]), 
       col = "red", pch = 16)


LM_Dist2Gene <- lm(log(Dist2Gene + 1) ~ as.factor(SelectInd), data = L1SingletonCoeffs)
summary(aov(LM_Dist2Gene))

LM_Freq <- lm(log(Freq) ~ as.factor(SelectInd), data = L1SingletonCoeffs)
summary(aov(LM_Freq))
LM_Freq_Sig <- lm(log(Freq) ~ as.factor(SelectInd), data = L1SingletonCoeffs,
                  subset = blnSig)
summary(aov(LM_Freq_Sig))

boxplot(log(Freq) ~ as.factor(SelectInd), data = L1SingletonCoeffs)

#######
# Regress against individual L1 positions
#######

L1SingletonCoeffs_subsetPenalized <- 
  subset(L1SingletonCoeffs, subset = robust.se > 0 & (!is.na(blnSelect)) &
           (!is.na(L1SingletonCoeffs$L1Start)) & (!is.na(L1SingletonCoeffs$L1End))&
           L1SingletonCoeffs$L1End >= 5900)
sum(L1SingletonCoeffs$L1End >= 5900, na.rm = T)
hist(L1SingletonCoeffs$L1End, seq(0, 6020, 10), xlim = c(5500, 6020))

# Create an indicator matrix for the presence of a L1 position inside an L1
L1End <- max(L1SingletonCoeffs$L1End, na.rm = T)
StWeights <- 1/L1SingletonCoeffs_subsetPenalized$robust.se / 
             mean(1/L1SingletonCoeffs_subsetPenalized$robust.se)
L1PosMat <- sapply(1:L1End, function(x){
  (L1SingletonCoeffs_subsetPenalized$L1Start >= x) & (L1SingletonCoeffs_subsetPenalized$L1End >= x)
})
dim(L1PosMat)

# Add an indicator for full-length at the end
L1PosMat <- cbind(L1PosMat, 1*L1SingletonCoeffs_subsetPenalized$blnFull)
L1PosMat <- 1*L1PosMat
L1PosMat[1:4, 1:5]
which(is.na(L1PosMat), arr.ind = T)

# Perform lasso and ridge regression
cat("Performing regularized regression ....")
GLM_Lasso <- glmnet(x = L1PosMat, y = 1*L1SingletonCoeffs_subsetPenalized$blnSelect,
                    alpha = 0.99, family = "binomial", weights = StWeights)
# GLM_Ridge <- glmnet(x = L1PosMat, y = 1*L1SingletonCoeffs_subsetPenalized$blnSelect,
#                     alpha = 0, family = "binomial", weights = StWeights)
# CV_Ridge <- cv.glmnet(x = L1PosMat, y = L1SingletonCoeffs_subsetPenalized$blnSelect,
#                alpha = 0, weights = StWeights)

# Coef_Lasso  <- coef(GLM_Lasso, CV_Lasso$lambda.1se)
# Coef_LassoV <- as.vector(Coef_Lasso)
# plot(Coef_LassoV)
# sum(Coef_LassoV > 0)
# Coef_Ridge <- coef(GLM_Ridge,
#    s = cv.glmnet(x = L1PosMat, y = L1SingletonCoeffs_subsetPenalized$blnSelect,
#                     alpha = 0)$lambda.1se)
# Predict_Ridge <- predict(GLM_Ridge, newx = L1PosMat,
#                    s = cv.glmnet(x = L1PosMat, y = L1SingletonCoeffs_subsetPenalized$blnSelect,
#                                  alpha = 0)$lambda.1se)
Predict_Lasso <- predict(GLM_Lasso, newx = L1PosMat,
                         s = cv.glmnet(x = L1PosMat, y = L1SingletonCoeffs_subsetPenalized$blnSelect,
                                       alpha = 0.99)$lambda.1se)

# Plot observed and predicted proportion of L1 
PropSelFull <- mean(L1SingletonCoeffs$blnSelect[L1SingletonCoeffs$blnFull], na.rm = T)
par(mfrow = c(1, 2), mar =  c(2, 2, 2, 1), oma = c(3, 3, 1, 0.1))
XLimLong  <- 6000
XLimShort <- 50
Cols <- rainbow(3)
plot(PropSmoothed_InsL$x, PropSmoothed_InsL$y, xlab = "",
     ylab = "", type = "l", xlim = c(0, XLimLong),
     ylim = c(0, 0.13), col = Cols[1], main = "A")
StartOrder <- order(L1SingletonCoeffs_subsetPenalized$L1Start)
#PredictP   <- exp(Predict_Ridge)/(1+exp(Predict_Ridge))
PredictP   <- exp(Predict_Lasso)/(1+exp(Predict_Lasso))
lines(L1SingletonCoeffs_subset$L1Start[InsLorder], 
      LM_All_Interact_binom$fitted.values[InsLorder],
      col = Cols[3])
lines(L1SingletonCoeffs_subsetPenalized$L1Start[StartOrder], PredictP[StartOrder],
      col = Cols[2])
arrows(XLimLong / 4, PropSelFull, 0, PropSelFull, length = 0.1)

plot(PropSmoothed_InsL$x, PropSmoothed_InsL$y, xlab = "",
      type = "l", xlim = c(0, XLimShort),
     ylim = c(0, 0.13), col = Cols[1], main = "B")
arrows(XLimShort / 4, PropSelFull, 0, PropSelFull, length = 0.1)
lines(L1SingletonCoeffs_subset$L1Start[InsLorder], 
      LM_All_Interact_binom$fitted.values[InsLorder],
      col = Cols[3])
lines(L1SingletonCoeffs_subsetPenalized$L1Start[StartOrder], PredictP[StartOrder],
      col = Cols[2])

mtext(text = "L1 insertion start [bp]", side = 1, outer = T, line = 1)
mtext(text = "Proportion of L1 with positive selection", side = 2, 
      outer = T, line = 1)
CreateDisplayPdf('D:/L1polymORF/Figures/PropSelectVsL1Start.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 4, width = 5)

# Coef_RidgeV <- as.vector(Coef_Ridge)
# sum(Coef_RidgeV > 0)
# sum(Coef_RidgeV > 0)
# plot(Coef_RidgeV, ylim = c(-0.002, 0.004), type = "l",
#      xlab = "Position of L1 present", ylab = "Effect on selection")
# lines(c(0, 10^4), c(0, 0), lty = 2, col = "red")

cat("done!\n")


# Check proportion selected among
L1SingletonCoeffs_with5P <- L1SingletonCoeffs[L1SingletonCoeffs$L1Start <= 10,]
table(L1SingletonCoeffs_with5P$blnFull, L1SingletonCoeffs_with5P$blnSelect)
fisher.test(L1SingletonCoeffs_with5P$blnFull, L1SingletonCoeffs_with5P$blnSelect)
sum(L1SingletonCoeffs_with5P$blnFull & L1SingletonCoeffs_with5P$blnSelect, na.rm = T)

##########################################
#                                        #
#       Plot frequency quantiles        #
#                                        #
##########################################

##########################
#                        #
#    Sample distances    #
#                        #
##########################

# Function to calculate p-values for distances to select features
L1Subset <- L1SingletonCoeffs
SelectDistP <- function(L1Subset, NrSample = 10^4){
  
  # Sample distances to various features of selected L1
  NrSelect <- sum(L1Subset$blnSelect)
  
  SampledMeandistMat <- sapply(1:NrSample, function(x){
    idx <- sample(1:nrow(L1Subset), NrSelect)
    c(
      Dist2Enhancer       = mean(L1Subset$Dist2Enhancer[idx]),
      Dist2Dist2EnhTxProm = mean(L1Subset$Dist2EnhTxProm[idx]),
      Dist2Txn = mean(L1Subset$Dist2Txn[idx]),
      Dist2Prom = mean(L1Subset$Dist2Prom[idx]),
      Dist2Repr = mean(L1Subset$Dist2Repr[idx]),
      LogDist2Repr = mean(log(L1Subset$Dist2Repr[idx] + 0.1)),
      Dist2Tfb = mean(L1Subset$Dist2Tfb[idx]),
      Dist2VistaOrGene = mean(L1Subset$Dist2VistaOrGene[idx])
    )
  })
  blnSel <- L1Subset$blnSelect
  MeanSelectDists <- c(
    Dist2Enhancer       = mean(L1Subset$Dist2Enhancer[blnSel]),
    Dist2Dist2EnhTxProm = mean(L1Subset$Dist2EnhTxProm[blnSel]),
    Dist2Txn = mean(L1Subset$Dist2Txn[blnSel]),
    Dist2Prom = mean(L1Subset$Dist2Prom[blnSel]),
    Dist2Repr = mean(L1Subset$Dist2Repr[blnSel]),
    LogDist2Repr = mean(log(L1Subset$Dist2Repr[blnSel] + 0.1)),
    Dist2Tfb = mean(L1Subset$Dist2Tfb[blnSel]),
    Dist2VistaOrGene = mean(L1Subset$Dist2VistaOrGene[blnSel])
  )
  rowSums(SampledMeandistMat <= MeanSelectDists) / NrSample
  
}
# Function to calculate p-values for distances to select features
FullDistP <- function(L1Subset, NrSample = 10^4){
  
  # Sample distances to various features of selected L1
  NrFull <- sum(L1Subset$blnFull)
  
  SampledMeandistMat <- sapply(1:NrSample, function(x){
    idx <- sample(1:nrow(L1Subset), NrFull)
    c(
      Dist2Enhancer       = mean(L1Subset$Dist2Enhancer[idx]),
      Dist2Dist2EnhTxProm = mean(L1Subset$Dist2EnhTxProm[idx]),
      Dist2Txn = mean(L1Subset$Dist2Txn[idx]),
      Dist2Prom = mean(L1Subset$Dist2Prom[idx]),
      Dist2Repr = mean(L1Subset$Dist2Repr[idx]),
      LogDist2Repr = mean(log(L1Subset$Dist2Repr[idx] + 0.1)),
      Dist2Tfb = mean(L1Subset$Dist2Tfb[idx])
    )
  })
  blnFull <- L1Subset$blnFull
  MeanFullDists <- c(
    Dist2Enhancer       = mean(L1Subset$Dist2Enhancer[blnFull]),
    Dist2Dist2EnhTxProm = mean(L1Subset$Dist2EnhTxProm[blnFull]),
    Dist2Txn = mean(L1Subset$Dist2Txn[blnFull]),
    Dist2Prom = mean(L1Subset$Dist2Prom[blnFull]),
    Dist2Repr = mean(L1Subset$Dist2Repr[blnFull]),
    LogDist2Repr = mean(log(L1Subset$Dist2Repr[blnFull] + 0.1)),
    Dist2Tfb = mean(L1Subset$Dist2Tfb[blnFull])
  )
  rowSums(SampledMeandistMat <= MeanFullDists) / NrSample
  
}


# Create a subset of coeffici
SDP_All   <- SelectDistP(L1SingletonCoeffs[!is.na(L1SingletonCoeffs$blnSelect),])
SDP_Fragm <- SelectDistP(L1SingletonCoeffs[which((!L1SingletonCoeffs$blnFull) &
                                                   !is.na(L1SingletonCoeffs$blnSelect)),])
SDP_Full  <- SelectDistP(L1SingletonCoeffs[which(L1SingletonCoeffs$blnFull),])
p.adjust(c(SDP_All, SDP_Full))
p.adjust(SDP_All)
p.adjust(SDP_Full)

FDP_All   <- FullDistP(L1SingletonCoeffs)


# Plot distance to closest regulatory feature for each combination of the 
# classification slected-non-selected, fragment-full
L1SingletonCoeffs$blnFull[is.na(L1SingletonCoeffs$blnFull)] <- FALSE
L1SingletonCoeffs$SelectFull <- paste(L1SingletonCoeffs$blnFull, L1SingletonCoeffs$blnSelect)
boxplot(log10(Dist2Repr) ~ SelectFull, data  = L1SingletonCoeffs)

L1SingletonCoeffs$Dist2Repr[L1SingletonCoeffs$blnFull & L1SingletonCoeffs$blnSelect]

# Calculate mean and standard deviation and plot
MeanDistBySelectFull <- aggregate(log10(L1SingletonCoeffs$Dist2Repr + 1), 
          by = list(L1SingletonCoeffs$blnSelect, 
                    L1SingletonCoeffs$blnFull),
          FUN = mean)
StErrDistBySelectFull <- aggregate(log10(L1SingletonCoeffs$Dist2Repr + 1), 
           by = list(L1SingletonCoeffs$blnSelect, L1SingletonCoeffs$blnFull),
           FUN = function(x) sqrt(var(x, na.rm = T) / length(x)))

bp <- barplot(MeanDistBySelectFull$x, ylim = c(0, 4.2), yaxt = "n",
              ylab = "Distance to repressed region [bp]")
Exps <- paste("c(", paste(paste("expression(10^", 0:4, ")"), collapse = ", "), ")")
axis(2, at = 0:4, eval(parse(text = Exps)))
axis(1, at = bp[,1], LETTERS[1:4])
AddErrorBars(MidX = bp[,1], MidY = MeanDistBySelectFull$x, TipWidth = 0.1,
             ErrorRange = StErrDistBySelectFull$x)
segments(x0 = bp[c(3, 3, 4),1], y0 = c(3.8, 4, 4),
         x1 = bp[c(3, 4, 4),1], y1 = c(4, 4, 3.8))
text(x = mean(bp[c(3, 4),1]), y = 4.15, "*", cex = 1.5)
CreateDisplayPdf('D:/L1polymORF/Figures/Dist2RegVsGroup.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 5, width = 5)

# Calculate mean distance to transcription factor binding site
MeanDist2TbBySelectFull <- aggregate(L1SingletonCoeffs$Dist2Tfb, 
                                     by = list(L1SingletonCoeffs$blnSelect, 
                                               L1SingletonCoeffs$blnFull),
                                     FUN = mean)

###################################
#                                 #
#    Fisher test for closeness    #
#                                 #
###################################

# Function to conduct fisher test for closeness
FisherTestClose <- function(blnSelect, Dist, CutOff = 2000){
  blnClose <- Dist <= CutOff
  fisher.test(blnSelect, blnClose)
}

# Fisher test for different distances
FisherTestClose(L1SingletonCoeffs$blnSelect, L1SingletonCoeffs$Dist2EnhancerOrGene)
FisherTestClose(L1SingletonCoeffs$blnSelect, L1SingletonCoeffs$Dist2Tfb, CutOff = 2000)

# Test for enrichment of different Tfbs
DistCutoff <- 1000
L1SingletonCoeffs_GRextended <- GRanges(seqnames = seqnames(L1SingletonCoeffs_GR),
   IRanges(start = start(L1SingletonCoeffs_GR) - DistCutoff,
           end = end(L1SingletonCoeffs_GR) + DistCutoff))
OL <- findOverlaps(L1SingletonCoeffs_GRextended, TfbGR)
idxSelect <- which(L1SingletonCoeffs$blnSelect)
blnSel <- OL@from %in% idxSelect

TfbNames <- unique(TfbdTable$name)

TfbEnrichP <- lapply(TfbNames, function(x){
  blnTfb <- TfbdTable$name[OL@to] == x
  if(sum(blnTfb & blnSel) > 1){
    fisher.test(blnTfb, blnSel)$p.value
  }
})
names(TfbEnrichP) <- TfbNames
TfbEnrichP <- unlist(TfbEnrichP)
min(p.adjust(TfbEnrichP))

###################################
#                                 #
#    Calculate population frequency        #
#                                     #
###################################


L1SingletonCoeffs$AFR
FreqPerPop <- aggregate.data.frame(
   L1SingletonCoeffs[,c("EUR", "EAS", "AMR", "SAS", "AFR")], 
   by = list(L1SingletonCoeffs$SelectInd), FUN = mean)

PSel_pos <- as.vector(FreqPerPop[1,-1] / FreqPerPop[2,-1])
PSel_neg <- as.vector(FreqPerPop[3,-1] / FreqPerPop[2,-1])
plot(t(PSel_pos), t(PSel_neg))
cor.test(t(PSel_pos), t(PSel_neg))$p.value
aggregate(Freq ~ SelectInd, data = L1SingletonCoeffs, FUN = mean)

CorPerPop <- sapply(as.character(UniquePops), function(x) {
  CT <- cor.test(L1SingletonCoeffs$SelectInd, L1SingletonCoeffs[,x],
           method = "spearman")
  c(Cor = CT$estimate, P = CT$p.value)
})
colnames(CorPerPop) <- UniquePops
CorPerPop
CT <- cor.test(L1SingletonCoeffs$coef, L1SingletonCoeffs[,"AFR"])
cor.test(L1SingletonCoeffs$coef, L1SingletonCoeffs[,"AFR"])

SelectPerChrom <- table(L1SingletonCoeffs$SelectInd, L1SingletonCoeffs$Chrom)
plot(SelectPerChrom[1,]/SelectPerChrom[2,],
     SelectPerChrom[3,]/SelectPerChrom[2,])

fisher.test(L1SingletonCoeffs$blnSig, L1SingletonCoeffs$Dist2Gene == 0)

###################################
#                                 #
#    Explore gene ontology        #
#                                 #
###################################

# Load the annotation package org.Hs.eg.db
library(org.Hs.eg.db)

# Create extended ranges of L1 insertions
DistCutoff <- 1000
L1SingletonCoeffs_GRextended <- 
  GRanges(seqnames = seqnames(L1SingletonCoeffs_GR),
          IRanges(start = start(L1SingletonCoeffs_GR) - DistCutoff,
                  end = end(L1SingletonCoeffs_GR) + DistCutoff))
# Get genomic ranges of genes
GeneGR <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
GeneGR_SelectPos  <- subsetByOverlaps(GeneGR, 
                       L1SingletonCoeffs_GRextended[which(L1SingletonCoeffs$blnSelect)])
GeneGR_SelectNeg  <- subsetByOverlaps(GeneGR, 
                                      L1SingletonCoeffs_GRextended[which(L1SingletonCoeffs$blnNotSelect)])
GeneGR_SelectAll  <- subsetByOverlaps(GeneGR, 
                                      L1SingletonCoeffs_GRextended[which(L1SingletonCoeffs$blnSig)])
GeneGR_NotSelect  <- subsetByOverlaps(GeneGR, 
                                      L1SingletonCoeffs_GRextended[which(!L1SingletonCoeffs$blnSig)])
GeneGR_All  <- subsetByOverlaps(GeneGR, 
                                      L1SingletonCoeffs_GRextended[which(!L1SingletonCoeffs$blnSig)])

# Write out gene ids for different genelists
writeLines(GeneGR_SelectPos@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_SelectPos")
writeLines(GeneGR_SelectNeg@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_SelectNeg")
writeLines(GeneGR_SelectAll@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_SelectAll")
writeLines(GeneGR_NotSelect@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_NotSelect")
writeLines(GeneGR_All@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_All")


# Get the gene symbol and name associated with the gene IDs in
# SelectGenesGR
GeneGROL@elementMetadata@listData$gene_id
cols <- c("SYMBOL", "GENENAME")
select(org.Hs.eg.db, keys = GeneGROL@elementMetadata@listData$gene_id,
       columns=cols, keytype="ENTREZID")

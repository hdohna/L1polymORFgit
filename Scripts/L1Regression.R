# The script below counts L1 and genome features per window and regresses L1 counts
# against the genome features

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
library(coin)

##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

# Specify file paths
DataPath           <- 'D:/L1polymORF/Data/'
G1000SamplePath    <- 'D:/L1polymORF/Data/1000GenomeSampleInfo.txt'
GapPath            <- 'D:/L1polymORF/Data/Gap_hg19.txt'
L1GRPath           <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
L1RefRangePath     <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
ChrLPath           <- 'D:/L1polymORF/Data/ChromLengthsHg19.Rdata'
GRangesSummaryPath <- 'D:/L1polymORF/Data/GRangesSummaries_1Mb.RData'

# Vector of predictor variables
# PredVars <- c("StemScores", "NotStemScores", "L1Count", "L1Count_Full",
#                "L1StartNum", "InsLength", "blnFull", "GC", "TargetFreq", 
#                "TeloFreq", "CpGCount", "length",  "cpgNum", "gcNum", 
#                "perCpg", "perGc", "obsExp", "decodeAvg", "decodeFemale", 
#                "decodeMale", "marshfieldAvg", "marshfieldMale", "genethonAvg",
#                "genethonFemale", "genethonMale", "phastCons", "EnhancerCount", 
#                "TxnCount", "PromCount", "ReprCount", "GeneCount") 
PredVars <- c("StemScores", "NotStemScores", "GC", "TargetFreq", 
                "TeloFreq", "CpGCount", "decodeAvg", "phastCons") 

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

# Load previously generated objects
load(L1GRPath)
load(ChrLPath)
load(L1RefRangePath)
load(GRangesSummaryPath)

cat("done\n")

########################################################
#                                                      #
#   Regress L1 count against data per genomic range    #
#                                                      #
########################################################

cat("Regress L1 counts ...")

# Create transparent point color
PCol  <- rgb(0,0,0, alpha = 0.1)

# Get a histogram of L1 count
hist(DataPerSummaryGR$L1Count)
plot(DataPerSummaryGR$L1Count - DataPerSummaryGR$L1Count_Full,
     DataPerSummaryGR$L1Count_Full, col = PCol, pch = 16)
cor.test(DataPerSummaryGR$L1Count - DataPerSummaryGR$L1Count_Full,
     DataPerSummaryGR$L1Count_Full)

# Regress LINE-1 count of 1000 genome against summaries
Form <- as.formula(paste("L1Count ~ ", paste(PredVars, collapse = "+")))
GLM_1000G <- glm(Form, data = DataPerSummaryGR, family = poisson)
summary(GLM_1000G)

# Regress LINE-1 count of reference genome against summaries
Form <- as.formula(paste("L1Count_ref ~ ", paste(PredVars, collapse = "+")))
GLM_Ref <- glm(Form, data = DataPerSummaryGR, family = poisson)
summary(GLM_Ref)

# Regress LINE-1 count of reference genome against summaries
GLM_L1_All <- glm(L1Count ~ NotStemScores + StemScores + GC +
                           CpGCount + EnhancerCount + TxnCount + 
                           + ReprCount + GeneCount + decodeAvg +
                           TargetFreq,
                         data = DataPerSummaryGR, family = poisson)
SUM_GLM_L1_All <- summary(GLM_L1_All)
SUM_GLM_L1_All$coefficients[,'Pr(>|z|)']
p.adjust(SUM_GLM_L1_All$coefficients[,'Pr(>|z|)'])
plot(DataPerSummaryGR$NotStemScores, DataPerSummaryGR$L1Count, col = PCol)
plot(DataPerSummaryGR$StemScores, DataPerSummaryGR$L1Count, col = PCol)

# Regress full-length L1 against summary variables
GLM_L1_Full <- glm(L1Count_Full ~ NotStemScores + StemScores + GC +
                     CpGCount + EnhancerCount + TxnCount + 
                     + ReprCount + GeneCount + decodeAvg +
                     TargetFreq,
                  data = DataPerSummaryGR, family = poisson)
SUM_GLM_L1_Full <- summary(GLM_L1_Full)
SUM_GLM_L1_Full
cbind(SUM_GLM_L1_Full$coefficients[,c('Estimate', 'Pr(>|z|)')], 
      SUM_GLM_L1_All$coefficients[,c('Estimate','Pr(>|z|)')])
plot(SUM_GLM_L1_Full$coefficients[,'Estimate'], SUM_GLM_L1_All$coefficients[,'Estimate'])

# Regress reference LINE-1 count against summaries
GLM_L1_Ref <- glm(L1Count_ref ~ NotStemScores + StemScores + GC +
                    CpGCount +  GeneCount + decodeAvg +
                    TxnCount + 
                    + ReprCount +TargetFreq + phastCons,
                  data = DataPerSummaryGR, family = poisson)
summary(GLM_L1_Ref)
SUM_GLM_L1_Ref <- summary(GLM_L1_Ref)
SUM_GLM_L1_Ref$coefficients[,'Pr(>|z|)']
p.adjust(SUM_GLM_L1_Ref$coefficients[,'Pr(>|z|)'])

# Compare coefficients for reference and non-reference L1
RefNonRefCoeff <- cbind(SUM_GLM_L1_Ref$coefficients[,'Estimate'], 
                        SUM_GLM_L1_All$coefficients[,'Estimate'])
RefNonRefCoeff <- RefNonRefCoeff / rowMeans(abs(RefNonRefCoeff))

plot(RefNonRefCoeff[-1, 1], RefNonRefCoeff[-1, 2])

# Check whether there are coefficients that are significant but opposite 
# directions for reference and non-reference L1
blnSigNonRef <- p.adjust(SUM_GLM_L1_All$coefficients[,'Pr(>|z|)']) < 0.1
blnSigNRef   <- p.adjust(SUM_GLM_L1_Ref$coefficients[,'Pr(>|z|)']) < 0.1
blnSigDiff <- blnSigNonRef & blnSigNRef & 
  (sign(SUM_GLM_L1_Ref$coefficients[,'Estimate']) != 
   sign(SUM_GLM_L1_All$coefficients[,'Estimate']) )
SUM_GLM_L1_Ref$coefficients[blnSigDiff,'Estimate']
SUM_GLM_L1_All$coefficients[blnSigDiff,'Estimate']


# Check GC content alone
glm(L1Count_ref ~ GC, data = DataPerSummaryGR, family = poisson)
glm(L1Count ~ GC, data = DataPerSummaryGR, family = poisson)
L1GRanges

cat("done\n")

########################################################
#                                                      #
#      Correlate L1 density with L1 frequency          #
#                                                      #
########################################################

cat("Correlate L1 density with L1 frequency ...")

# Get L1 frequency per summary range
OL <- findOverlaps(SummaryGR, L1_1000G_GR_hg19)
MeanFreq <- aggregate(L1_1000G$Frequency[OL@to], by = list(OL@from), FUN = mean)
DataPerSummaryGR$MeanFreq <- NA
DataPerSummaryGR$MeanFreq[MeanFreq$Group.1] <- MeanFreq$x

# Get predicted L1 frequency
L1Predict <- predict(GLM_L1_All, newdata = DataPerSummaryGR)

# Test for correlation
cor.test(DataPerSummaryGR$MeanFreq, exp(L1Predict))
cor.test(DataPerSummaryGR$MeanFreq, L1Predict, method = "spearman")
spearman_test(DataPerSummaryGR$MeanFreq ~ L1Predict)
spearman_test(DataPerSummaryGR$MeanFreq ~ DataPerSummaryGR$L1Count)

plot(DataPerSummaryGR$MeanFreq, DataPerSummaryGR$L1Count, col = PCol, pch = 16)
plot(DataPerSummaryGR$L1Count, DataPerSummaryGR$MeanFreq, col = PCol, pch = 16)
MeanFreqPerCount <- aggregate(MeanFreq ~ L1Count, data = DataPerSummaryGR,
                              FUN = mean)
plot(MeanFreqPerCount$L1Count, MeanFreqPerCount$MeanFreq)
VarFreqPerCount <- aggregate(MeanFreq ~ L1Count, data = DataPerSummaryGR,
                              FUN = var)
boxplot(MeanFreq ~ L1Count, data = DataPerSummaryGR)

cat("done\n")

#################################################################
#                                                               #
#   Different regression for high medium and low frequency L1   #
#                                                               #
#################################################################

cat("Regress L1 counts for different L1 density classes...")

# Regress low frequency LINE-1 count against summaries
GLM_L1_low <- glm(L1Count_low ~ NotStemScores + StemScores + GC +
                    CpGCount + EnhancerCount + TxnCount + 
                    + ReprCount + GeneCount + decodeAvg +
                    TargetFreq,
                  data = DataPerSummaryGR, family = poisson)
SUM_GLM_L1_low <- summary(GLM_L1_low)

# Regress medium frequency LINE-1 count against summaries
GLM_L1_med <- glm(L1Count_med ~ NotStemScores + StemScores + GC +
                    CpGCount + EnhancerCount + TxnCount + 
                    + ReprCount + GeneCount + decodeAvg +
                    TargetFreq,
                  data = DataPerSummaryGR, family = poisson)
SUM_GLM_L1_med <- summary(GLM_L1_med)

# Regress medium frequency LINE-1 count against summaries
GLM_L1_high <- glm(L1Count_high ~ NotStemScores + StemScores + GC +
                    CpGCount + EnhancerCount + TxnCount + 
                    + ReprCount + GeneCount + decodeAvg +
                    TargetFreq,
                  data = DataPerSummaryGR, family = poisson)
SUM_GLM_L1_high <- summary(GLM_L1_high)


# Put regression coefficients for low, medium and high density L1 into a matrix
CoeffMat <- cbind(SUM_GLM_L1_low$coefficients[,c('Estimate')], 
      SUM_GLM_L1_med$coefficients[,c('Estimate')],
      SUM_GLM_L1_high$coefficients[,c('Estimate')])
CoeffMat <- CoeffMat - rowMeans(CoeffMat)
CoeffMat <- CoeffMat / sqrt(apply(CoeffMat, 1, var))
colnames(CoeffMat) <- c("LowFreq", "MedFreq", "HighFreq")

# Plot lines for all coefficients (except intercept)
Cols <- rainbow(nrow(CoeffMat) - 1)
plot(c(1, 3), c(min(CoeffMat), max(CoeffMat)), type = "n")
for (i in 2:nrow(CoeffMat)){
  lines(CoeffMat[i,], col = Cols[i - 1])
}
L1FreqClass <- 1:3
FreqChange <- sapply(2:nrow(CoeffMat), function(i){
  LM <- lm(CoeffMat[i,] ~ L1FreqClass)
  summary(LM)$coefficients[2,c('Estimate', 'Pr(>|t|)')]
})
colnames(FreqChange) <- rownames(CoeffMat)[-1]

cat("done\n")

#################################################################
#                                                               #
#              GEE regression                 #
#                                                               #
#################################################################

# cat("Performing GEE regression ...")

# Regress LINE-1 count against DNAse scores
GLM_L1_DNAse_gee <- gee(L1Count ~ NotStemScores + StemScores + GC +
                          CpGCount + EnhancerCount + TxnCount + 
                          + ReprCount + GeneCount + decodeAvg +
                          TargetFreq,
                         data = DataPerSummaryGR, family = "poisson",
                         corstr =  "exchangeable", Mv = 1,
                        id = ChrNum)
SUM_GLM_L1_DNAse_gee <- summary(GLM_L1_DNAse_gee)
coef(SUM_GLM_L1_DNAse_gee)[,'Estimate']
se <- SUM_GLM_L1_DNAse_gee$coefficients[, "Robust S.E."]
Ps <- 1 - pnorm(abs(coef(SUM_GLM_L1_DNAse_gee)[,'Estimate']) / se)
format(Ps, digits = 22)
cbind(coef(SUM_GLM_L1_DNAse_gee)[,'Estimate'] - se * qnorm(0.9995),
      coef(SUM_GLM_L1_DNAse_gee)[,'Estimate'] + se * qnorm(0.9995))

GLM_L1_DNAse_stem <- glm(L1Count ~ StemScores, data = DataPerSummaryGR,
                    family = "poisson")
summary(GLM_L1_DNAse_stem)
GLM_L1_DNAse_NotStem <- glm(L1Count ~ NotStemScores, data = DataPerSummaryGR,
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

# ##############################################
# #                                            #
# #   Analyze difference between               #
# #   expected and observed heterozygosity     #
# #                                            #
# ##############################################
# 
# # Get expected and observed number of homozygotes
# L1_1000G$HomoExp <- length(SampleColumns) * L1_1000G$Frequency^2
# L1_1000G$HomoObs <- rowSums(L1_1000G[,SampleColumns] == 2)
# L1_1000G$HomoDiff <- (L1_1000G$HomoObs - L1_1000G$HomoExp) 
# mean(L1_1000G$HomoDiff, na.rm = T)
# hist(L1_1000G$HomoDiff, breaks = -500:200)
# hist(L1_1000G$HomoExp, breaks = 0:2000, xlim = c(0, 10)) 
# 
# # Check whether differemce between expected and observed homozygosity differes
# # between L1 in genes and outside genes
# boxplot(HomoDiff ~ blnOLGene, data = L1_1000G)
# t.test(HomoDiff ~ blnOLGene, data = L1_1000G)
# wilcox.test(HomoDiff ~ blnOLGene, data = L1_1000G)
# with(L1_1000G, mean(HomoDiff[blnOLGene], na.rm = T))
# with(L1_1000G, mean(HomoDiff[!blnOLGene], na.rm = T))
# 
# ##########################################
# #                                        #
# #   Test gene-singleton association      #
# #                                        #
# ##########################################
# 
# # Create boolean variable for intersection with genes
# L1SingletonCoeffs$blnOLGene <- L1SingletonCoeffs$Dist2Gene == 0
# fisher.test(L1SingletonCoeffs$blnOLGene, L1SingletonCoeffs$blnNotSelect)
# wilcox.test(CoefSt ~ blnOLGene, data = L1SingletonCoeffs)
# t.test(CoefSt ~ blnOLGene, data = L1SingletonCoeffs)
# aggregate(CoefSt ~ blnOLGene, data = L1SingletonCoeffs, FUN = mean)
# 
# wilcox.test(coef ~ blnOLGene, data = L1SingletonCoeffs)
# t.test(coef ~ blnOLGene, data = L1SingletonCoeffs)
# aggregate(coef ~ blnOLGene, data = L1SingletonCoeffs, FUN = mean)
# 
# sum(L1_1000G$blnOLGene) / 3060

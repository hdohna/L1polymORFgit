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
RecombDataPath  <- "D:/L1polymORF/Data/hg19RecombRate.txt"

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

# Read in table with regulatory elements
RegTable      <- read.table("D:/L1polymORF/Data/EncodeBroadHMM/ChromHMMcombined.txt",
                            header = T)
blnAllCellTypes <- RegTable$CellType == 
  "Gm12878,H1hesc,Hepg2,Hmec,Hsmm,Huvec,K562,Nhek,Nhlf"

# Get indices of different regulatory elements
idxEnhancer <- grep("Enhancer", RegTable$name)
idxTxn      <- grep("Txn", RegTable$name)
idxProm     <- grep("Promoter", RegTable$name)
idxHetero   <- grep("Heterochrom", RegTable$name)
idxRepr     <- union(grep("Insulator", RegTable$name),
                     grep("Repressed", RegTable$name))

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
EnhancerGR  <- makeGRangesFromDataFrame(RegTable[idxEnhancer,])
TxnGR       <- makeGRangesFromDataFrame(RegTable[idxTxn,])
PromGR      <- makeGRangesFromDataFrame(RegTable[idxProm, ])
ReprGR      <- makeGRangesFromDataFrame(RegTable[idxRepr, ])

# Read in table with L1HS from the referemce genome
L1RefTab <- read.csv(L1RefPath)
L1RefGR <- makeGRangesFromDataFrame(L1RefTab, seqnames.field = "genoName",
                                    start.field = "genoStart",
                                    end.field = "genoEnd")

# Read in CpG data and turn it into GRanges
CpGtable <- read.delim(CpGPath, header = T)
CpGGR    <- makeGRangesFromDataFrame(CpGtable, start.field = "chromStart",
                                     end.field = "chromEnd")

# Read in file and create GRanges
RecData <- read.delim("D:/L1polymORF/Data/hg19RecombRate.txt")
Rec_GR  <- makeGRangesFromDataFrame(RecData)

# Add columns to 1000 genome data
L1_1000G$blnOLGene  <- overlapsAny(L1_1000G_GR_hg19, GeneGR)
L1_1000G$blnOLProm  <- overlapsAny(L1_1000G_GR_hg19, PromGR)
L1_1000G$blnOLExon  <- overlapsAny(L1_1000G_GR_hg19, ExonGR)
L1_1000G$blnOLCpG   <- overlapsAny(L1_1000G_GR_hg19, CpGGR)
L1_1000G$Dist2CpG   <- Dist2Closest(L1_1000G_GR_hg19, CpGGR)
L1_1000G$L1StartNum <- as.numeric(as.character(L1_1000G$L1Start))
L1_1000G$L1EndNum   <- as.numeric(as.character(L1_1000G$L1End))
L1_1000G$blnFull    <- L1_1000G$L1StartNum <= 1 & L1_1000G$L1EndNum >= 6000

cat("done!\n")


###################################################
#                                                 #
#     Summarize data per genomic range            #
#                                                 #
###################################################

cat("Summarizing data per genomic range ...")
# Summarize DNAse hypesensitivity data per RangeWidth
DNAse2Summarize  <- cbind(DNAseData[,c("chrom", "start", "end")], 
                          ScoreByStem)
DNAseSummaryList <- AggregateDataPerGRanges(
  Data2Summarize = DNAse2Summarize,
  Cols2Summarize = colnames(ScoreByStem),
  SeqNameCol = "chrom",
  StartCol = "start",
  EndCol = "end",
  RangeWidth = RangeWidth, 
  ChromLengths = ChromLengthsHg19,
  blnAddGRInfo = T,
  blnReturnDFOnly = F,
  NoValue = 0
)

# Create summary genomic ranges and data per summary range
SummaryGR        <- DNAseSummaryList$SummaryGR
DataPerSummaryGR <- DNAseSummaryList$SummarizedData
rm(list = "DNAseSummaryList")
gc()

###################
# Add more variables to summary
###################

# Add overlap counts
DataPerSummaryGR$L1Count <- countOverlaps(SummaryGR, L1_1000G_GR_hg19)
DataPerSummaryGR$L1Count_Full <- countOverlaps(SummaryGR, 
                                  L1_1000G_GR_hg19[which(L1_1000G$blnFull)])

# Add other L1 info to summary
L1Data2Summarize <- L1_1000G[ ,c("CHROM", "POS", "L1StartNum", "InsLength",
                                 "blnFull")]
L1Data2Summarize$CHROM <- paste("chr", L1Data2Summarize$CHROM, sep = "")
L1SummaryDf <- AggregateDataPerGRanges(
  Data2Summarize = L1Data2Summarize,
  SummaryGR = SummaryGR,
  Cols2Summarize = c("L1StartNum", "InsLength", "blnFull"),
  SeqNameCol = "CHROM",
  StartCol = "POS",
  EndCol = "POS",
  RangeWidth = RangeWidth, 
  ChromLengths = ChromLengthsHg19,
  blnAddGRInfo = T)

# Merge existing data with new ones
DataPerSummaryGR <- merge(DataPerSummaryGR, L1SummaryDf, sort = F)
rm(list = "L1SummaryDf")
gc()


# Add GC to summary
cat("calculating GC content ...")
SummaryViews <- BSgenomeViews(BSgenome.Hsapiens.UCSC.hg19, SummaryGR)
GCFreq       <- letterFrequency(SummaryViews, letters = "GC")
DataPerSummaryGR$GC <- GCFreq
cat("done\n")

# Add target site frequency to summary
cat("calculating target site frequency ...")
SeqWindows <- getSeq(BSgenome.Hsapiens.UCSC.hg19, SummaryGR)
TTTTAAFreq   <- vcountPattern("TTTTAA", SeqWindows)
TTAAAAFreq   <- vcountPattern("TTAAAA", SeqWindows)
DataPerSummaryGR$TargetFreq <- TTTTAAFreq + TTAAAAFreq
cat("done\n")

# Add CpG data
DataPerSummaryGR$CpGCount <- countOverlaps(SummaryGR, CpGGR)
CpGSummaryDF <- AggregateDataPerGRanges(
  Data2Summarize = CpGtable,
  Cols2Summarize = c("length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp"),
  SeqNameCol     = "chrom",
  StartCol       = "chromStart",
  EndCol         = "chromEnd",
  RangeWidth     = RangeWidth, 
  ChromLengths = ChromLengthsHg19,
  blnAddGRInfo = T)
CpGSummaryDF[is.na(CpGSummaryDF$length),c("length", "cpgNum", "gcNum")] <- 0

# Merge existing data with new ones
DataPerSummaryGR <- merge(DataPerSummaryGR, CpGSummaryDF, sort = F)
rm(list = "CpGSummaryDF")
gc()

# Add recombination data
RecSummaryDF <- AggregateDataPerGRanges(
  Data2Summarize = RecData,
  Cols2Summarize = c("decodeAvg",  "decodeFemale", "decodeMale", 
                     "marshfieldAvg", "marshfieldMale", "genethonAvg",
                     "genethonFemale", "genethonMale"),
  SeqNameCol     = "chrom",
  StartCol       = "chromStart",
  EndCol         = "chromEnd",
  RangeWidth     = RangeWidth, 
  ChromLengths = ChromLengthsHg19,
  blnAddGRInfo = T)

# Merge existing data with new ones
DataPerSummaryGR <- merge(DataPerSummaryGR, RecSummaryDF, sort = F)
rm(list = "CpGSummaryDF")
gc()

# Add various counts (better: calculate sum of width)
CoverGR <- GeneGR
CalcGRCoverage <- function(SummaryGR, CoverGR){
  WidthDF <- data.frame(Chrom = as.vector(seqnames(CoverGR)),
    Start = start(CoverGR), End = end(CoverGR), Width = width(CoverGR))
  AggWidth <- AggregateDataPerGRanges(
    SummaryGR = SummaryGR,
    Data2Summarize = WidthDF,
    Cols2Summarize = c("Width"),
    SeqNameCol     = "Chrom",
    StartCol       = "Start",
    EndCol         = "End",
    ChromLengths = ChromLengthsHg19,
    SummaryFun = sum,
    blnAddGRInfo = T)
  AggWidth$Width
}
DataPerSummaryGR$EnhancerCount <- CalcGRCoverage(SummaryGR, EnhancerGR)
DataPerSummaryGR$TxnCount      <- CalcGRCoverage(SummaryGR, TxnGR)
DataPerSummaryGR$PromCount     <- CalcGRCoverage(SummaryGR, PromGR)
DataPerSummaryGR$ReprCount     <- CalcGRCoverage(SummaryGR, ReprGR)
DataPerSummaryGR$GeneCount     <- CalcGRCoverage(SummaryGR, GeneGR)

cat("done!\n")

########################################################
#                                                      #
#   Regress L1 count against data per genomic range    #
#                                                      #
########################################################

# Create transparent point color
PCol  <- rgb(0,0,0, alpha = 0.1)

# Get a histogram of L1 count
hist(DataPerSummaryGR$L1Count)
plot(DataPerSummaryGR$L1Count - DataPerSummaryGR$L1Count_Full,
     DataPerSummaryGR$L1Count_Full, col = PCol, pch = 16)
cor.test(DataPerSummaryGR$L1Count - DataPerSummaryGR$L1Count_Full,
     DataPerSummaryGR$L1Count_Full)

# Regress LINE-1 count against summaries
GLM_L1_All <- glm(L1Count ~ NotStemScores + StemScores + GC +
                           CpGCount + EnhancerCount + TxnCount + 
                           + ReprCount + GeneCount + decodeAvg +
                    TargetFreq,
                         data = DataPerSummaryGR, family = poisson)
SUM_GLM_L1_All <- summary(GLM_L1_All)
SUM_GLM_L1_All$coefficients[,'Pr(>|z|)']
exp(SUM_GLM_L1_DNAse_both$coefficients[,'Estimate'])
plot(DataPerSummaryGR$NotStemScores, DNAseSummary$L1Count, col = PCol)
plot(DataPerSummaryGR$NotStemScores, DNAseSummary$L1Count, col = PCol)

GLM_L1_Full <- glm(L1Count_Full ~ NotStemScores + StemScores + GC +
                    CpGCount + EnhancerCount + TxnCount + 
                    + ReprCount + GeneCount +  decodeAvg,
                  data = DataPerSummaryGR, family = poisson)
SUM_GLM_L1_Full <- summary(GLM_L1_Full)
SUM_GLM_L1_Full
cbind(SUM_GLM_L1_Full$coefficients[,c('Estimate', 'Pr(>|z|)')], 
      SUM_GLM_L1_All$coefficients[,c('Estimate','Pr(>|z|)')])
plot(SUM_GLM_L1_Full$coefficients[,'Estimate'], SUM_GLM_L1_All$coefficients[,'Estimate'])

# # Regress LINE-1 count against DNAse scores
# GLM_L1_DNAse_gee <- gee(L1Count ~ NotStemScores + StemScores + GC, 
#                          data = DataPerSummaryGR, family = "poisson",
#                          corstr =  "exchangeable", Mv = 1, 
#                         id = ChrNum)
# SUM_GLM_L1_DNAse_gee <- summary(GLM_L1_DNAse_gee)
# coef(SUM_GLM_L1_DNAse_gee)[,'Estimate']
# se <- SUM_GLM_L1_DNAse_gee$coefficients[, "Robust S.E."]
# Ps <- 1 - pnorm(abs(coef(SUM_GLM_L1_DNAse_gee)[,'Estimate']) / se)
# format(Ps, digits = 22)
# cbind(coef(SUM_GLM_L1_DNAse_gee)[,'Estimate'] - se * qnorm(0.9995),
#       coef(SUM_GLM_L1_DNAse_gee)[,'Estimate'] + se * qnorm(0.9995))
# 
# GLM_L1_DNAse_stem <- glm(L1Count ~ StemScores, data = DataPerSummaryGR,
#                     family = "poisson")
# summary(GLM_L1_DNAse_stem)
# GLM_L1_DNAse_NotStem <- glm(L1Count ~ NotStemScores, data = DataPerSummaryGR,
#                     family = "poisson")
# GLM_L1_DNAse_GC <- glm(L1Count ~ GC, data = DNAseSummary,
#                             family = "poisson")
# summary(GLM_L1_DNAse_GC)
# GLM_L1_DNAse_fragm <- glm(L1Count_fragm ~ StemScores + NotStemScores, data = DNAseSummary,
#                     family = "poisson")
# summary(GLM_L1_DNAse_fragm)
# GLM_L1_DNAse_full <- glm(L1Count_full ~ StemScores + NotStemScores, data = DNAseSummary,
#                           family = "poisson")
# summary(GLM_L1_DNAse_full)
# max(DNAseSummary$L1Count_full)
# 
# plot(DNAseSummary$StemScores, DNAseSummary$NotStemScores, col = PCol)
# plot(DNAseSummary$StemScores, DNAseSummary$L1Count, col = PCol)
# plot(DNAseSummary$NotStemScores, DNAseSummary$L1Count, col = PCol)
# 
# cor(DNAseSummary$StemScores, DNAseSummary$NotStemScores)
# cor(DNAseSummary$StemScores, DNAseSummary$GC)
# cor(DNAseSummary$NotStemScores, DNAseSummary$GC)
# 
# 
# # Determine proportion of L1 overlapping with heterochromatin
# # blnOLHetero <- overlapsAny(L1_1000G_GR_hg19, HeteroGR)
# # mean(blnOLHetero)
# # sum(width(HeteroGR)/10^6) / sum(ChromLengthsHg19/10^6)
# # 
# # fisher.test(blnOLHetero, L1_1000G$InsLength >= 6000)
# 
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

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

# Number of info columns in vcf file
NrInfoCols   <- 9

# False discovery rate for selected L1
FDR <- 0.1

# Specify range width for DNAse analysis
RangeWidth <- 10^5

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
PromGR <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene)

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

##########################################
#                                        #
#   Regress intersection with genes      #
#                                        #
##########################################

#######
# Regress against L1 start
#######

# Indicator variable for intersection with genes
L1_1000G$blnOLGene  <- overlapsAny(L1_1000G_GR_hg19, GeneGR)
L1_1000G$blnOLProm  <- overlapsAny(L1_1000G_GR_hg19, PromGR)
L1_1000G$L1StartNum <- as.numeric(as.character(L1_1000G$L1Start))
L1_1000G$L1EndNum   <- as.numeric(as.character(L1_1000G$L1End))
L1_1000G$blnFull    <- L1_1000G$L1StartNum <= 1 & L1_1000G$L1EndNum >= 6000
hist(L1_1000G$L1StartNum, breaks = 0:6100)
sum(L1_1000G$L1StartNum == 2, na.rm = T)
sum(L1_1000G$blnOLProm)

# Create transparent point color
PCol  <- rgb(0,0,0, alpha = 0.1)
PCol1 <- rgb(1,0,0, alpha = 0.1)
PCol2 <- rgb(0,0,1, alpha = 0.1)

# Regress indicator for gene overlap against L1 start abd frequency
L1_1000G_subset <- L1_1000G[L1_1000G$Frequency >= 0.01, ]
L1_1000G_subset <- L1_1000G
LogRegGeneOLVsL1Start <- glm(blnOLGene ~ L1StartNum + Frequency + blnFull, 
                             family = "binomial", data = L1_1000G_subset)
summary(LogRegGeneOLVsL1Start)

table(L1_1000G$L1StartNum)[1:20]

# Regress indicator for promoter overlap against L1 start abd frequency
LogRegPromOLVsL1Start <- glm(blnOLProm ~ L1StartNum  + Frequency, 
                             family = "binomial", data = L1_1000G)
summary(LogRegPromOLVsL1Start)

# Plot proportion gene overlap against L1 start
L1OLVsL1StartSmoothed <- supsmu(L1_1000G_subset$L1StartNum,  1*L1_1000G_subset$blnOLGene)
plot(L1OLVsL1StartSmoothed$x, L1OLVsL1StartSmoothed$y, type = "l", xlab = "L1 start",
     ylab = "Proportion of L1s in genes")
points(L1_1000G$L1StartNum[L1_1000G$blnOLGene], 
       rep(0.31, sum(L1_1000G$blnOLGene)), col = PCol1)
points(L1_1000G$L1StartNum[!L1_1000G$blnOLGene], 
       rep(0.32, sum(!L1_1000G$blnOLGene)), col = PCol2)
plot(L1OLVsL1StartSmoothed$x, L1OLVsL1StartSmoothed$y, type = "l", xlab = "L1 start",
     ylab = "Proportion of L1s in genes", xlim = c(0, 30),
     ylim = c(0.38, 0.55))
points(L1_1000G$L1StartNum[L1_1000G$blnOLGene], 
       rep(0.42, sum(L1_1000G$blnOLGene)), col = PCol1, pch = 15, cex = 1.3)
points(L1_1000G$L1StartNum[!L1_1000G$blnOLGene], 
       rep(0.41, sum(!L1_1000G$blnOLGene)), col = PCol2, pch = 15, cex = 1.3)

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

##########################################
#                                        #
#   Regress strand alignedness           #
#                                        #
##########################################

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

# Plot proportion gene overlap against L1 start
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
length(GeneGR_All)

# Check whether gene size determines overlap
GeneOLCount <- countOverlaps(GeneGR, L1_1000G_GR_hg19)
boxplot(log(width(GeneGR)) ~ GeneOLCount)

# Plot smoothed proportion overlap against gene length
OLVsGeneLSmoothed <- supsmu(width(GeneGR),  GeneOLCount)
OLVsGeneLSmoothed_log10 <- supsmu(log10(width(GeneGR)),  GeneOLCount)
plot(log10(width(GeneGR)),  GeneOLCount, col = PCol, pch = 16,
     xlab = "Gene length")
lines(OLVsGeneLSmoothed$x, OLVsGeneLSmoothed$y, col = "red")
plot(width(GeneGR), GeneOLCount, xlab = "Gene length", xlim = c(0, 2*10^6),
     col = PCol, pch = 16)
lines(OLVsGeneLSmoothed$x, OLVsGeneLSmoothed$y, col = "red")

sort(width(GeneGR), decreasing = T)[1:5]
max(1*GeneOLCount)

# Write out gene ids for different genelists
writeLines(GeneGR_All@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_All")
writeLines(GeneGR_LowFreq@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_LowFreq")
writeLines(GeneGR_HiFreq@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_HiFreq")

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

###################################################
#                                                 #
#         Analyze enrichment of GO terms          #
#                                                 #
###################################################

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
                        idxGeneGR = idxUniProt)

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
                   PDZ + Concanavalin + Ligand + MAM + PTP + PDEase + HostReceptor, 
                   data = GeneTable,
                   family = "poisson")
Sum_GLM_OLCount <- summary(GLM_OLCount)
Padj <- p.adjust(Sum_GLM_OLCount$coefficients[,'Pr(>|z|)'])
cbind(Sum_GLM_OLCount$coefficients[Padj < 0.05, 'Estimate'], Padj[Padj < 0.05])

# Create 

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
L1_1000G$blnOLGlyco  <- overlapsAny(L1_1000G_GR_hg19, GeneGR[idxGlyco])
LogRegGlycoOLVsL1Start <- glm(blnOLGlyco ~ L1StartNum + Frequency + blnFull, 
                             family = "binomial", data = L1_1000G)
summary(LogRegGlycoOLVsL1Start)
L1GlycoOLVsL1StartSmoothed <- supsmu(L1_1000G$L1StartNum,  1*L1_1000G$blnOLGlyco)
plot(L1GlycoOLVsL1StartSmoothed$x, L1GlycoOLVsL1StartSmoothed$y, type = "l", 
     xlab = "L1 start", ylab = "Proportion of L1s in glycoproteins")
L1_1000G_noGlyco <- L1_1000G[!L1_1000G$blnOLGlyco, ]
L1noGlycoOLVsL1StartSmoothed <- supsmu(L1_1000G_noGlyco$L1StartNum,  1*L1_1000G_noGlyco$blnOLGene)
plot(L1noGlycoOLVsL1StartSmoothed$x, L1noGlycoOLVsL1StartSmoothed$y, type = "l", 
     xlab = "L1 start", ylab = "Proportion of L1s in genes")

# Get entrez IDs of glycoproteins 
GeneGR[idxGlyco]
keytypes(org.Hs.eg.db)
select(org.Hs.eg.db,
       keys = GeneGR@elementMetadata@listData$gene_id[idxGlycoL1],
       columns = c("SYMBOL", "GENENAME"), keytype = "ENTREZID")

###################################################
#                                                 #
#     Predict L1 insertion based on DNAse         #
#                                                 #
###################################################

# Function to create genomic ranges for a particular chromosome number
CreateGR <- function(i){
  CL <- ChromLengthsHg19[i]
  GRStarts <- seq(1, CL, RangeWidth)
  GREnds  <- c(GRStarts[-1] - 1, CL)
  GRanges(seqnames = names(CL), IRanges(start = GRStarts, end = GREnds))
}

# Create genomic ranges for 
SummaryGR <- CreateGR(1)
for (i in 2:length(ChromLengthsHg19)){
  NewGR <- CreateGR(i)
  SummaryGR <- c(SummaryGR, NewGR)
}

# Summarize score by stem-cells and non-stem cells (move to separate script)
# ScoreByStem <- t(sapply(1:nrow(DNAseData), function(i){
#   IDs <- strsplit(as.character(DNAseData$sourceIds[i]), ",")[[1]]
#   Scores <- as.numeric(strsplit(as.character(DNAseData$sourceScores[i]), ",")[[1]])
#   c(StemScores = sum(Scores[IDs %in% StemIDs]),
#     NotStemScores = sum(Scores[IDs %in% NotStemIDs]))
# }))

# Create genomic ranges for the DNAse data and then aggregate by SummaryGR
DNAseGR <- makeGRangesFromDataFrame(DNAseData)
DNAseSummaryOL <- findOverlaps(DNAseGR, SummaryGR)

# Calculate mean DNAse per summary genomic range
DNAseBySummaryGR <- aggregate.data.frame(ScoreByStem[DNAseSummaryOL@from, ], 
                                         by = list(DNAseSummaryOL@to), FUN = mean)
DNAseSummary <- as.data.frame(matrix(0, nrow = length(SummaryGR), ncol = 2))
DNAseSummary[DNAseBySummaryGR$Group.1, ] <- DNAseBySummaryGR[,2:3]
colnames(DNAseSummary) <- colnames(DNAseBySummaryGR)[2:3]
DNAseSummary$L1Count   <- countOverlaps(SummaryGR, L1_1000G_GR_hg19)
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
GLM_L1_DNAse_sum <- glm(L1Count ~ ScoreSum + pStem, data = DNAseSummary,
                         family = "poisson")
summary(GLM_L1_DNAse_sum)
GLM_L1_DNAse_both <- glm(L1Count ~ NotStemScores + StemScores, data = DNAseSummary,
                    family = "poisson")
summary(GLM_L1_DNAse_both)
GLM_L1_DNAse_stem <- glm(L1Count ~ StemScores, data = DNAseSummary,
                    family = "poisson")
summary(GLM_L1_DNAse_stem)
GLM_L1_DNAse_NotStem <- glm(L1Count ~ NotStemScores, data = DNAseSummary,
                    family = "poisson")
summary(GLM_L1_DNAse_NotStem)
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


# Determine proportion of L1 overlapping with heterochromatin
# blnOLHetero <- overlapsAny(L1_1000G_GR_hg19, HeteroGR)
# mean(blnOLHetero)
# sum(width(HeteroGR)/10^6) / sum(ChromLengthsHg19/10^6)
# 
# fisher.test(blnOLHetero, L1_1000G$InsLength >= 6000)

    
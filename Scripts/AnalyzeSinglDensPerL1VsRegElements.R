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

# Read in table with Info about 1000 genome samples 
SampleInfo_1000Genome <- read.table(G1000SamplePath, as.is = T, header = T)
table(SampleInfo_1000Genome$super_pop)


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
idxEnhancer    <- grep("Enhancer", RegTable$name)
idxTxn         <- grep("Txn", RegTable$name)
idxProm        <- grep("Promoter", RegTable$name)
idxHetero      <- grep("Heterochrom", RegTable$name)
idxRepr        <- union(grep("Insulator", RegTable$name),
                        grep("Repressed", RegTable$name))

NT <- table(RegTable$name)
blnAllCellTypes <- RegTable$CellType == 
                     "Gm12878,H1hesc,Hepg2,Hmec,Hsmm,Huvec,K562,Nhek,Nhlf"
RegTableAll <- RegTable[blnAllCellTypes,]
sort(table(RegTable$CellType))

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
RegGR    <- makeGRangesFromDataFrame(RegTable)
RegGRAll <- makeGRangesFromDataFrame(RegTableAll)
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

cat("done!\n")

##########################################
#                                        #
#        Add columns                     #
#                                        #
##########################################

# Indicator for full-length
L1SingletonCoeffs$blnFull <- L1SingletonCoeffs$InsLength >= 6000
sum(L1SingletonCoeffs$InsLength <= 100)

# Indicator for significant effect
L1SingletonCoeffs$blnSig <- p.adjust(L1SingletonCoeffs$Pr...z..) < 0.05

# Indicator for  selection
L1SingletonCoeffs$blnSelect <- L1SingletonCoeffs$blnSig &
  L1SingletonCoeffs$coef < 0
L1SingletonCoeffs$blnSelect[is.na(L1SingletonCoeffs$blnSelect)] <- FALSE
  
table(L1SingletonCoeffs$blnSelect, L1SingletonCoeffs$blnFull)
fisher.test(L1SingletonCoeffs$blnSelect, L1SingletonCoeffs$blnFull)
nrow(L1SingletonCoeffs)

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

##########################################
#                                        #
#        Regress coefficients            #
#                                        #
##########################################


#######
# Regress against length
#######

# Subset L1 coefficients
sum(is.na(L1SingletonCoeffs$InsLength))
L1SingletonCoeffs_subset <- subset(L1SingletonCoeffs, 
                                   subset = robust.se > 0 & (!is.na(blnSelect)))
LM_All_Interact_binom <- glm(1L*blnSelect ~ InsLength + blnFull, 
                             data = L1SingletonCoeffs_subset, 
                             weights = 1/robust.se,
                             family = quasibinomial)
summary(LM_All_Interact_binom)
mean(1/L1SingletonCoeffs_subset$robust.se)

# Smooth proportions of 
PropSmoothed_InsL <- supsmu(L1SingletonCoeffs_subset$InsLength,
                            1*L1SingletonCoeffs_subset$blnSelect,
                            wt = 1/L1SingletonCoeffs_subset$robust.se)
plot(PropSmoothed_InsL$x, PropSmoothed_InsL$y, xlab = "L1 insertion length [bp]",
     ylab = "Proportion of L1 with positive selection signal", type = "l",
     ylim = c(0, 0.11))
lines(PropSmoothed_InsL$x[PropSmoothed_InsL$x >=6000], 
      PropSmoothed_InsL$y[PropSmoothed_InsL$x>=6000])
InsLorder <- order(L1SingletonCoeffs_subset$InsLength)
lines(L1SingletonCoeffs_subset$InsLength[InsLorder], 
      LM_All_Interact_binom$fitted.values[InsLorder], lty = 2)
CreateDisplayPdf('D:/L1polymORF/Figures/PropSelectVsInsLength.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 5, width = 5)

##########################
#                        #
#    Sample distances    #
#                        #
##########################

# Function to calculate p-values for distances to select features
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
      Dist2Tfb = mean(L1Subset$Dist2Tfb[idx])
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
    Dist2Tfb = mean(L1Subset$Dist2Tfb[blnSel])
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
SDP_All   <- SelectDistP(L1SingletonCoeffs)
SDP_Fragm <- SelectDistP(L1SingletonCoeffs[which(!L1SingletonCoeffs$blnFull),])
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
DistCutoff <- 2000
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

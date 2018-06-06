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

# Read in table with known L1 
L1Catalogue <- read.csv("D:/L1polymORF/Data/L1CatalogExtended.csv", as.is = T)
L1Catalogue$Allele[is.na(L1Catalogue$Allele)] <- 1

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
blnInRef1         <- (L1Catalogue$end_HG38 - L1Catalogue$start_HG38) > 6000 
L1CatalogL1Mapped <- L1Catalogue[blnL1Mapped & blnAllele1,]

# Get strand from 1000 genome data
Strand1000G <- sapply(as.character(L1_1000G[,8]), GetStrandFrom1000GVcf)

# Create genomic ranges for catalog L1
L1CatalogGR <- GRanges(seqnames = L1CatalogL1Mapped$Chromosome,
                       ranges = IRanges(start = pmin(L1CatalogL1Mapped$start_HG38,
                                                     L1CatalogL1Mapped$end_HG38),
                                        end = pmax(L1CatalogL1Mapped$start_HG38,
                                                   L1CatalogL1Mapped$end_HG38)),
                       strand = L1CatalogL1Mapped$strand_L1toRef)

DistObj        <- distanceToNearest(L1CatalogGR, L1_1000G_GRList_hg38$LiftedRanges, 
                             ignore.strand = T) 
DistCat2_1000G  <- DistObj@elementMetadata@listData$distance

# Get indices of 1000 Genome and catalog elements that match
idx1000G <- nearest(L1CatalogGR, L1_1000G_GRList_hg38$LiftedRanges)
L1_1000G_GRList_hg38$idxUniqueMapped
L1CatalogMatch1000G <- L1CatalogL1Mapped[DistObj@from[DistCat2_1000G < 100], ]
idx1000GMatchCat    <- L1_1000G_GRList_hg38$idxUniqueMapped[DistObj@to[DistCat2_1000G < 100]]
L1CatalogGRMatched  <- L1CatalogGR[DistObj@from[DistCat2_1000G < 100]]
sum(DistCat2_1000G < 100)
sum(DistCat2_1000G < 10)
L1_1000G$InsLength[idx1000GMatchCat]
L1_1000G$InsLength[idx1000GMatchCat]

blnSameStrand <- L1CatalogMatch1000G$strand_L1toRef == Strand1000G[idx1000GMatchCat]
boxplot(L1_1000G$InsLength[idx1000GMatchCat] ~ blnSameStrand)
t.test(L1_1000G$InsLength[idx1000GMatchCat] ~ blnSameStrand)
table(L1CatalogMatch1000G$strand_L1toRef, Strand1000G[idx1000GMatchCat])

# add a column with dummy activity sums
L1CatalogMatch1000G$ActivityDummy <- 1

# Get rows of L1_1000G table that can be matched to catalog
L1_1000G_match <- L1_1000G[idx1000GMatchCat, ]
L1_1000G_match$InsLength
L1CatalogMatch1000G$Activity[L1_1000G_match$InsLength < 500]

# Get expected and observed activity sums
ExpectedAct <- 2 * L1CatalogMatch1000G$ActivityNum %*% L1_1000G_match$Frequency
ObservedAct <- sapply(SampleColumns, function(x){
  L1CatalogMatch1000G$ActivityNum %*% L1_1000G_match[,x]})
hist(ObservedAct)
mean(ObservedAct)
length(ObservedAct)
mean(L1_1000G_reduced$Frequency)

############
#  Process histone data
############

# Read table with open compartment data
HistoTableFiles <- list.files("D:/L1polymORF/Data/", 
                              pattern = "H3K", full.names = T)

##########
#  Process ChromHMM
##########

# Read in table with regulatory elements
RegTable       <- read.table("D:/L1polymORF/Data/ChromHMM", header = T)
# RegTable      <- read.table("D:/L1polymORF/Data/EncodeBroadHMM/wgEncodeBroadHmmH1hescHMM.txt",
#                              col.names = c("bin", "chrom", "ChromStart",
#                                            "ChromEnd", "name", "score", "strand",  
#                                            "thickStart", "thickEnd", "itemRgb"))


idxEnhancer    <- grep("Enhancer", RegTable$name)
idxTxnTable    <- grep("Txn", RegTable$name)
idxHeteroTable <- grep("Heterochrom", RegTable$name)

NT <- table(RegTable$name)

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
RegGR <- makeGRangesFromDataFrame(RegTable, start.field = "ChromStart", 
                                  end.field = "ChromEnd")
EnhancerGR  <- RegGR[idxEnhancer]
TxnGR       <- RegGR[idxTxnTable]
RegNonHetGR <- RegGR[-idxHeteroTable]

##########
#  Process ChromHMM
##########

TfbdTable <- read.table("D:/L1polymORF/Data/wgEncodeRegTfbsClusteredV3.txt",
               col.names = c("bin", "chrom", "start", "end", "name", "score",
                             "expCount", "expNums", "expScores"))
TfbGR <- makeGRangesFromDataFrame(TfbdTable)

##########
# Process loops
##########

# Auxiliary unction to create genomic ranges for both sides of a loop
getLoopGRs <- function(FileName){
  Loops <- read.delim(FileName)
  Loops$chr1 <- paste("chr", Loops$chr1, sep = "")
  Loops$chr2 <- paste("chr", Loops$chr2, sep = "")
  LoopsGR1 <- makeGRangesFromDataFrame(Loops, seqnames.field = "chr1", 
                                       start.field="x1", end.field = "x2")
  LoopsGR2 <- makeGRangesFromDataFrame(Loops, seqnames.field = "chr2", 
                                       start.field="y1", end.field = "y2")
  blnOverlapLoop <- overlapsAny(LoopsGR1, LoopsGR2)
  AllLoops <- c(LoopsGR1, LoopsGR2[!blnOverlapLoop])
  GRangesList(LoopsGR1 = LoopsGR1, LoopsGR2 = LoopsGR2, AllLoops = AllLoops)
}

getLoopDomainGRs <- function(FileName){
  Loops <- read.delim(FileName)
  Loops$chr1 <- paste("chr", Loops$chr1, sep = "")
  Loops$chr2 <- paste("chr", Loops$chr2, sep = "")
  makeGRangesFromDataFrame(Loops, seqnames.field = "chr1", 
                           start.field="x1", end.field = "y2")
}

# List of files with domains and loops
DomainFiles <- list.files(HiCFolderPath, pattern = "domainlist.txt",
                          full.names = T)
DomainFiles <- DomainFiles[! DomainFiles %in% grep(".gz", DomainFiles, value = T)]
LoopFiles <- list.files(HiCFolderPath, pattern = "looplist.txt",
                        full.names = T)
LoopFiles <- LoopFiles[! LoopFiles %in% grep(".gz", LoopFiles, value = T)]

# Get list of domain and loop ranges
DomainGRList   <- lapply(DomainFiles, getLoopGRs)
LoopGRList     <- lapply(LoopFiles, getLoopGRs)
LoopDomainGRList     <- lapply(LoopFiles, getLoopDomainGRs)
names(DomainGRList) <- sapply(DomainFiles, 
                              function(x) strsplit(x, "_")[[1]][2])
names(LoopGRList) <- sapply(LoopFiles, 
                            function(x) strsplit(x, "_")[[1]][2])
names(LoopDomainGRList) <- sapply(LoopFiles, 
                                  function(x) strsplit(x, "_")[[1]][2])

# Summarize domain ranges (union or intersect)
DomainGR_Intersect <- DomainGRList[[1]]$AllLoops
LoopGR_Intersect   <- LoopGRList[[1]]$AllLoops
LoopDomainGR_Intersect   <- LoopDomainGRList[[1]]
for (i in 2:length(LoopGRList)){
  DomainGR_Intersect     <- intersect(DomainGR_Intersect, DomainGRList[[i]]$AllLoops)
  LoopGR_Intersect       <- intersect(LoopGR_Intersect, LoopGRList[[i]]$AllLoops)
  LoopDomainGR_Intersect <- intersect(LoopDomainGR_Intersect, LoopDomainGRList[[i]])
}
blnOverlap <- overlapsAny(L1SingletonCoeffs_GR, DomainGR_Intersect)
sum(blnOverlap)
table(blnOverlap, L1SingletonCoeffs$blnSelect)
table(blnOverlap[which(L1SingletonCoeffs$blnFull)], 
      L1SingletonCoeffs$blnSelect[which(L1SingletonCoeffs$blnFull)])

# Make genomic ranges for L1SingletonCoeffs
L1SingletonCoeffs$chromosome <- paste("chr", L1SingletonCoeffs$Chrom, sep = "")
L1SingletonCoeffs_GR <- makeGRangesFromDataFrame(L1SingletonCoeffs, 
                                                 seqnames.field = "chromosome",
                                                 start.field = "Pos",
                                                 end.field = "Pos")
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

# Calculate distance to closest enhancer
L1SingletonCoeffs$Dist2Tfb <- Dist2Closest(L1SingletonCoeffs_GR, TfbGR)

L1SingletonCoeffs$Dist2TfbTn <- pmin(L1SingletonCoeffs$Dist2Txn,
                                     L1SingletonCoeffs$Dist2Tfb)

# Distance to closest regulatory element
L1SingletonCoeffs$Dist2Reg <- Dist2Closest(L1SingletonCoeffs_GR, RegGR)
RegOverlap <- findOverlaps(L1SingletonCoeffs_GR, RegGR)
table(RegTable$name[RegOverlap@to])

# Distance to closest regulatory element
L1SingletonCoeffs$Dist2RegNonHet <- Dist2Closest(L1SingletonCoeffs_GR, RegNonHetGR)

# Closest regulatory element
idxNearestReg <- nearest(L1SingletonCoeffs_GR, RegNonHetGR)
RegNames <- RegTable$name[-idxHeteroTable]
L1SingletonCoeffs$ClosestReg <- RegNames[idxNearestReg]

# Determine distance to loops and domains
L1SingletonCoeffs$Dist2Loop <- Dist2Closest(L1SingletonCoeffs_GR, 
                                            LoopGR_Intersect)
L1SingletonCoeffs$Dist2LoopDomain <- Dist2Closest(L1SingletonCoeffs_GR, 
                                                  LoopDomainGR_Intersect)
L1SingletonCoeffs$FullDist2Loop <- L1SingletonCoeffs$Dist2Loop * L1SingletonCoeffs$blnFull
# L1SingletonCoeffs$Dist2Domain <- Dist2Closest(L1SingletonCoeffs_GR, 
#                                               DomainGR_Intersect)
L1SingletonCoeffs$Dist2LoopBins <- cut(L1SingletonCoeffs$Dist2Loop,
                                       breaks = seq(0, 24*10^6, 10^5))


HADf <- AggregateHistoneMarks(HistoTableFiles[1], L1SingletonCoeffs_GR) 
  
table(L1SingletonCoeffs$blnSelect, L1SingletonCoeffs$InsLBins)

##########################################
#                                        #
#        Plot coefficients               #
#                                        #
##########################################

# Plot coefficients vs insertion length
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
lines(PropSmoothed_InsL$x[PropSmoothed_InsL$x>=6000], 
      PropSmoothed_InsL$y[PropSmoothed_InsL$x>=6000])
InsLorder <- order(L1SingletonCoeffs_subset$InsLength)
lines(L1SingletonCoeffs_subset$InsLength[InsLorder], 
      LM_All_Interact_binom$fitted.values[InsLorder], lty = 2)
CreateDisplayPdf('D:/L1polymORF/Figures/PropSelectVsInsLength.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 5, width = 5)

#######
# Regress against distance to genes
#######

# Regress coefficient against distance to genes
LM_Full <- lm(coef ~ Dist2Gene, data = L1SingletonCoeffs, 
         subset = L1SingletonCoeffs$InsLength >= 6000 & L1SingletonCoeffs$robust.se > 0,
         weights = 1/L1SingletonCoeffs$robust.se)
summary(LM_Full)

LM_Fragm <- lm(coef ~ Dist2Gene, data = L1SingletonCoeffs, 
              subset = L1SingletonCoeffs$InsLength < 5900 & L1SingletonCoeffs$robust.se > 0,
              weights = 1/L1SingletonCoeffs$robust.se)
summary(LM_Fragm)

LM_All <- lm(coef ~ Dist2Gene, data = L1SingletonCoeffs, 
              subset = L1SingletonCoeffs$robust.se > 0,
              weights = 1/L1SingletonCoeffs$robust.se)
LM_All <- lm(coef ~ Dist2Gene, data = L1SingletonCoeffs)
summary(LM_All)
sum(L1SingletonCoeffs$robust.se == 0)

LM_All_Interact_binom <- glm(1*blnSelect ~ Dist2Gene + blnFull + blnFull*Dist2Gene, 
                             data = L1SingletonCoeffs, 
                             subset =robust.se > 0 & (!is.na(blnSelect)),
                             weights = 1/L1SingletonCoeffs$robust.se, family  = quasibinomial)
summary(LM_All_Interact_binom)



L1SingletonCoeffs[which(L1SingletonCoeffs$blnSelect),]

max(L1SingletonCoeffs$Dist2Loop)

#######
# Regress against distance to loops
#######

LM_InteractDist2Loop_binom <- glm(1L*blnSelect ~ InsLength + blnFull + Dist2Loop*InsLength, 
                             data = L1SingletonCoeffs_subset, 
                             weights = 1/robust.se,
                             family = quasibinomial)
summary(LM_InteractDist2Loop_binom)



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
    c(Dist2Gene = mean(L1Subset$Dist2Gene[idx]),
      Dist2Loop = mean(L1Subset$Dist2Loop[idx]),
      Dist2Enhancer = mean(L1Subset$Dist2Enhancer[idx]),
      Dist2EnhancerOrGene = mean(L1Subset$Dist2EnhancerOrGene[idx]),
      Dist2Tfb = mean(L1Subset$Dist2Tfb[idx]),
      Dist2Txn = mean(L1Subset$Dist2Txn[idx]),
      Dist2TfbTn = mean(L1Subset$Dist2TfbTn[idx]),
      Dist2Reg = mean(L1Subset$Dist2Reg[idx]),
      Dist2RegNonHet = mean(L1Subset$Dist2RegNonHet[idx]),
      LogDist2RegNonHet = mean(log(L1Subset$Dist2RegNonHet[idx] + 1))
    )
  })
  blnSel <- L1Subset$blnSelect
  MeanSelectDists <- c(Dist2Gene = mean(L1Subset$Dist2Gene[blnSel]),
                       Dist2Loop = mean(L1Subset$Dist2Loop[blnSel]),
                       Dist2Enhancer = mean(L1Subset$Dist2Enhancer[blnSel]),
                       Dist2EnhancerOrGene = mean(L1Subset$Dist2EnhancerOrGene[blnSel]),
                       Dist2Tfb = mean(L1Subset$Dist2Tfb[blnSel]),
                       Dist2Txn = mean(L1Subset$Dist2Txn[blnSel]),
                       Dist2TfbTn = mean(L1Subset$Dist2TfbTn[blnSel]),
                       Dist2Reg = mean(L1Subset$Dist2Reg[blnSel]),
                       Dist2RegNonHet = mean(L1Subset$Dist2RegNonHet[blnSel]),
                       LogDist2RegNonHet = mean(log(L1Subset$Dist2RegNonHet[blnSel] + 0.1))
                       
  )
  rowSums(SampledMeandistMat <= MeanSelectDists) / NrSample
  
}


table(L1SingletonCoeffs$Dist2Reg > 0, L1SingletonCoeffs$blnSelect)
table(L1SingletonCoeffs$ClosestReg[which(L1SingletonCoeffs$blnFull &
                                     L1SingletonCoeffs$blnSelect)])
table(L1SingletonCoeffs$ClosestReg[which(!L1SingletonCoeffs$blnFull &
                                           L1SingletonCoeffs$blnSelect)]) /
table(L1SingletonCoeffs$ClosestReg[which(!L1SingletonCoeffs$blnSelect)]) * 
  sum(!L1SingletonCoeffs$blnSelect) / sum(L1SingletonCoeffs$blnSelect) 
table(L1SingletonCoeffs$ClosestReg)
table(RegTable$name)

# Create a subset of coeffici
SDP_all <- SelectDistP(L1SingletonCoeffs)
p.adjust(SDP_all[c(1, 3, 4, 5, 6, 8)])
SelectDistP(L1SingletonCoeffs[which(!L1SingletonCoeffs$blnFull),])
SelectDistP(L1SingletonCoeffs[which(L1SingletonCoeffs$blnFull),])

# Plot distance to closest regulatory feature for each combination of the 
# classification slected-non-selected, fragment-full
L1SingletonCoeffs$blnFull[is.na(L1SingletonCoeffs$blnFull)] <- FALSE
L1SingletonCoeffs$SelectFull <- paste(L1SingletonCoeffs$blnFull, L1SingletonCoeffs$blnSelect)
boxplot(log10(Dist2RegNonHet) ~ SelectFull, data  = L1SingletonCoeffs)

# Calculate mean and standard deviation and plot
MeanDistBySelectFull <- aggregate(log10(L1SingletonCoeffs$Dist2RegNonHet + 1), 
          by = list(L1SingletonCoeffs$blnSelect, 
                    L1SingletonCoeffs$blnFull),
          FUN = mean)
StErrDistBySelectFull <- aggregate(log10(L1SingletonCoeffs$Dist2RegNonHet + 1), 
           by = list(L1SingletonCoeffs$blnSelect, L1SingletonCoeffs$blnFull),
           FUN = function(x) sqrt(var(x, na.rm = T) / length(x)))

bp <- barplot(MeanDistBySelectFull$x, ylim = c(0, 4.2), yaxt = "n",
              ylab = "Distance to regulatory element [bp]")
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

MeanDist2TbBySelectFull <- aggregate(L1SingletonCoeffs$Dist2Tfb, 
                                  by = list(L1SingletonCoeffs$blnSelect, 
                                            L1SingletonCoeffs$blnFull),
                                  FUN = mean)




t.test(L1SingletonCoeffs$Dist2Tfb [L1SingletonCoeffs$blnFull] ~ 
         L1SingletonCoeffs$blnSelect[L1SingletonCoeffs$blnFull])
L1SingletonCoeffs$Dist2Txn[L1SingletonCoeffs$blnSelect &
                             L1SingletonCoeffs$blnFull]
hist(log10(L1SingletonCoeffs$Dist2Txn[L1SingletonCoeffs$blnFull]))
(sum(L1SingletonCoeffs$Dist2Enhancer[L1SingletonCoeffs$blnFull] <= 40000, na.rm = T)/
  sum(L1SingletonCoeffs$blnFull, na.rm = T))^3
sum(L1SingletonCoeffs$blnFull & L1SingletonCoeffs$blnSelect, na.rm = T)
SampledMeanDist <- sapply(1:10000, function(x){
  idx <- sample(which(L1SingletonCoeffs$blnFull), 3)
  mean(L1SingletonCoeffs$Dist2Enhancer[idx])
})
sum(SampledMeanDist <= 
      mean(L1SingletonCoeffs$Dist2Enhancer[which(L1SingletonCoeffs$blnSelect &
             L1SingletonCoeffs$blnFull)])) / 10000

Distances  <- L1SingletonCoeffs$Dist2Txn
blnSelect  <- L1SingletonCoeffs$blnSelect
blnFull    <- L1SingletonCoeffs$blnFull
SampleSize = 10^5

#AnalyzeDistanceInteraction <- function(Distances, blnSelect, 
#   blnFull, SampleSize = 10000){
sum(is.na(L1SingletonCoeffs$Dist2Enhancer))
blnNA           <- is.na(Distances) | is.na(blnSelect) | is.na(blnFull)
Distances       <- Distances[!blnNA]
blnSelect       <- blnSelect[!blnNA]
blnFull         <- blnFull[!blnNA]

blnSelectFull   <- blnSelect[blnFull]
blnSelectFragm  <- blnSelect[!blnFull]
Distances_Full  <- Distances[blnFull]
Distances_Fragm <- Distances[!blnFull]

MeanDiffFull <- mean(Distances_Full[blnSelectFull]) - 
  mean(Distances_Full[!blnSelectFull])
MeanDiffFragm <- mean(Distances_Fragm[blnSelectFragm]) - 
  mean(Distances_Fragm[!blnSelectFragm])
MeanDiffBoth  <- MeanDiffFull - MeanDiffFragm
NrFull        <- sum(blnFull)
NrFragm       <- sum(!blnFull)
NrSelectFull  <- sum(blnSelectFull)
NrSelectFragm <- sum(blnSelectFragm)


SampledDiffMat <- sapply(1:SampleSize, function(x) {
  SampleSelectFull   <- sample(1:NrFull,  size = NrSelectFull)
  SampleSelectFragm  <- sample(1:NrFragm, size = NrSelectFragm)
  MeanDiffFull <- mean(Distances_Full[SampleSelectFull]) - 
    mean(Distances_Full[-SampleSelectFull])
  MeanDiffFragm <- mean(Distances_Fragm[SampleSelectFragm]) - 
    mean(Distances_Fragm[-SampleSelectFragm])
  MeanDiffBoth = MeanDiffFull - MeanDiffFragm
  c(MeanDiffFull = MeanDiffFull, MeanDiffFragm = MeanDiffFragm,
    MeanDiffBoth = MeanDiffBoth)
  
})
dim(SampledDiffMat)
sum(SampledDiffMat["MeanDiffFull", ]  <= MeanDiffFull) / SampleSize
sum(SampledDiffMat["MeanDiffFragm", ] <= MeanDiffFragm) / SampleSize
sum(SampledDiffMat["MeanDiffBoth", ]  <= MeanDiffBoth) / SampleSize

fisher.test(rbind(c(NrSelectFull, NrFull - NrSelectFull),
                  c(NrSelectFragm, NrFragm - NrSelectFragm)))

####################################################
#                                                  #
#   Correlate activity with selection coefficient  #
#                                                  #
####################################################

# Create vector of chromosome and positions to match L1s from different files
ChrPos1 <- paste(L1_1000G_match$CHROM, L1_1000G_match$POS)
ChrPos2 <- paste(L1SingletonCoeffs$Chrom, L1SingletonCoeffs$Pos)
L1SinglCoeffs_match <- L1SingletonCoeffs[match(ChrPos1, ChrPos2),]
L1SinglCoeffs_match$Activity <- L1CatalogMatch1000G$ActivityNum
plot(L1SinglCoeffs_match$coef, L1CatalogMatch1000G$ActivityNum)
LM_Activity_binom <- glm(blnSelect ~ Activity, 
                             data = L1SinglCoeffs_match, 
                             weights = 1/robust.se,
                             family = quasibinomial)
summary(LM_Activity_binom)
sum(L1SinglCoeffs_match$blnSelect)

pSelect <- mean(L1SingletonCoeffs$blnSelect[L1SingletonCoeffs$blnFull],
     na.rm = T)
(1 - pSelect)^nrow(L1SinglCoeffs_match)

##############################################
#                                            #
#    Test for L1 insertion length effects    #
#                                            #
##############################################

# Kernel smooth 
# blnNotNA <- !is.na(L1SingletonCoeffs$InsLength)
# KSM_NotSelect <- ksmooth(x = L1SingletonCoeffs$InsLength[which(blnNotNA &
#                                (!L1SingletonCoeffs$blnSelect))], 
#         y = L1SingletonCoeffs$Dist2Gene[which(blnNotNA &
#               (!L1SingletonCoeffs$blnSelect))], 
#         kernel = "normal", n.points = 25, bandwidth = 200)
# KSM_Select <- ksmooth(x = L1SingletonCoeffs$InsLength[which(blnNotNA &
#                             (L1SingletonCoeffs$blnSelect))], 
#                       y = L1SingletonCoeffs$Dist2Gene[which(blnNotNA &
#                               (L1SingletonCoeffs$blnSelect))], 
#                      kernel = "normal", n.points = 25, bandwidth = 200)
# sum(is.na(KSM_NotSelect$y))
# sum(is.na(KSM_Select$y))
# plot(KSM_Select$x, KSM_Select$y, type = "l", col = "red")
# lines(c(0, 10000), rep(mean(KSM_Select$y), 2), lty = 2, col = "red")
# lines(KSM_NotSelect$x, KSM_NotSelect$y)
# lines(c(0, 10000), rep(mean(KSM_NotSelect$y), 2), lty = 2)
# 
# plot(KSM_NotSelect$y - mean(KSM_NotSelect$y), KSM_Select$y - KSM_NotSelect$y)
# cor(KSM_NotSelect$y - mean(KSM_NotSelect$y), KSM_Select$y - KSM_NotSelect$y)

SmDistNonSelect <- supsmu(L1SingletonCoeffs$InsLength[!L1SingletonCoeffs$blnSelect], 
                          L1SingletonCoeffs$Dist2Gene[!L1SingletonCoeffs$blnSelect])
plot(L1SingletonCoeffs$InsLength, L1SingletonCoeffs$Dist2Enhancer)
points(L1SingletonCoeffs$InsLength[L1SingletonCoeffs$blnSelect], 
       L1SingletonCoeffs$Dist2Enhancer[L1SingletonCoeffs$blnSelect],
       pch = 16, col = "red")
lines(SmDistNonSelect$x, SmDistNonSelect$y, col = "blue", lwd = 2)
plot(L1SingletonCoeffs$Dist2Gene, L1SingletonCoeffs$Dist2Enhancer)
points(L1SingletonCoeffs$Dist2Gene[L1SingletonCoeffs$blnSelect], 
       L1SingletonCoeffs$Dist2Enhancer[L1SingletonCoeffs$blnSelect],
       pch = 16, col = "red")
cor.test(L1SingletonCoeffs$Dist2Gene, L1SingletonCoeffs$Dist2Enhancer)
SampledCors <- sapply(1:1000, function(x){
  SampleRows <- sample(1:nrow(L1SingletonCoeffs), size = sum(L1SingletonCoeffs$blnSelect,
                                                             na.rm = T))
  cor(L1SingletonCoeffs$Dist2Gene[SampleRows], 
      L1SingletonCoeffs$Dist2Enhancer[SampleRows])
  
})
sum(cor(L1SingletonCoeffs$Dist2Gene[which(L1SingletonCoeffs$blnSelect)], 
         L1SingletonCoeffs$Dist2Enhancer[which(L1SingletonCoeffs$blnSelect)]) >=
      SampledCors) / length(SampledCors)  
table(L1SingletonCoeffs$Dist2Enhancer <= 100, L1SingletonCoeffs$blnFull)
table(L1SingletonCoeffs$Dist2Enhancer <= 1000, L1SingletonCoeffs$blnSelect,
      L1SingletonCoeffs$blnFull)

t.test(L1SingletonCoeffs$Dist2EnhancerOrGene ~ L1SingletonCoeffs$blnSelect)
t.test(L1SingletonCoeffs$Dist2Gene ~ L1SingletonCoeffs$blnSelect)
t.test(L1SingletonCoeffs$Dist2Enhancer ~ L1SingletonCoeffs$blnSelect)

# Function to get mean, upper and lower quantile per insertion length class
QMat <- function(SampleMat, LowerQ = 0.025, UpperQ = 0.975,
                 Plot = T, xVals = NULL, Method = c("Mean", "Median")){
  QMat <- matrix(nrow = 3, ncol = nrow(SampleMat))
  QMat[1,]   <- rowMeans(SampleMat)
  #  QMat[1,]   <- rowMedians(SampleMat)
  QMat[2:3,] <- apply(SampleMat, 1, 
                      FUN = function(x) quantile(x, c(LowerQ, UpperQ)))
  if(Plot){
    idxR <- ncol(QMat):1
    polygon(c(xVals, xVals[idxR]), c(QMat[2, ], QMat[3, idxR]), 
            col = "grey", border = NA)
    
  }
  QMat
}


# Function to analyze selection patterns per L1 insertion length
SelectPatternPerInsL <- function(ResponseVarName = "Dist2Gene", 
  InsLBins = cut(L1SingletonCoeffs$InsLength, breaks = seq(0, 6500, 500)), 
  Method = c("Mean", "Median"),
  NrSamples = 1000){
  
  # Specify the number of samples and initialize matrices for sampled
  # quantities
  LengthBins <- unique(InsLBins)
  LengthBins <- LengthBins[!is.na(LengthBins)]
  SampledVar <- matrix(nrow = length(LengthBins), ncol = NrSamples)
  SampledCor <- rep(NA, NrSamples)
  
  if (Method[1] == "Mean"){ 
    MeanVal   <- mean(L1SingletonCoeffs[,ResponseVarName])
    # Loop to sample response variable
    for (i in 1:NrSamples){
      SampledBins <- sample(InsLBins, size = nrow(L1SingletonCoeffs))
      MeanPerL    <- aggregate(L1SingletonCoeffs[, ResponseVarName],
                               by = list(SampledBins), 
                               FUN = function(x) mean(x, na.rm = T))
      # Aggregate response variable by insertion length and selection signal
      MeanPerLSelect <- aggregate(L1SingletonCoeffs[, ResponseVarName],
                                  by = list(SampledBins, L1SingletonCoeffs$blnSelect), 
                                  FUN = function(x) mean(x, na.rm = T))
      # Difference between non-selected L1 and overall mean
      DiffNonSelectMean <- MeanPerLSelect$x[!MeanPerLSelect$Group.2] - 
        MeanVal
      
      # Difference between selected and non-selected L1 and overall mean
      SelectMatch  <- match(MeanPerLSelect$Group.1[MeanPerLSelect$Group.2],
                            MeanPerLSelect$Group.1[!MeanPerLSelect$Group.2])
      DiffSelectNonSelect <- MeanPerLSelect$x[MeanPerLSelect$Group.2] - 
        MeanPerLSelect$x[!MeanPerLSelect$Group.2][SelectMatch]
      SampledVar[,i] <- MeanPerL$x
      SampledCor[i] <- cor(DiffNonSelectMean[SelectMatch], DiffSelectNonSelect)
    }
    # Aggregate response variable by insertion length and selection signal
    MeanPerL <- aggregate(L1SingletonCoeffs[,c("InsLength", ResponseVarName)],
                              by = list(InsLBins), 
                              FUN = function(x) mean(x, na.rm = T))
    # Aggregate response variable by insertion length and selection signal
    MeanPerLSelect <- aggregate(L1SingletonCoeffs[,c("InsLength", ResponseVarName)],
                            by = list(InsLBins, L1SingletonCoeffs$blnSelect), 
                            FUN = function(x) mean(x, na.rm = T))
    
  } else {
    MeanVal   <- median(L1SingletonCoeffs[,ResponseVarName])
    for (i in 1:NrSamples){
      SampledBins <- sample(InsLBins, size = nrow(L1SingletonCoeffs))
      MeanPerL    <- aggregate(L1SingletonCoeffs[,ResponseVarName],
                               by = list(SampledBins), 
                               FUN = function(x) median(x, na.rm = T))
      SampledVar[,i] <- MeanPerL$x
    }
    # Aggregate response variable by insertion length and selection signal
    MeanPerL <- aggregate(L1SingletonCoeffs[,c("InsLength", ResponseVarName)],
                          by = list(InsLBins), 
                          FUN = function(x) median(x, na.rm = T))
    # Aggregate response variable by insertion length and selection signal
    MeanPerLSelect <- aggregate(L1SingletonCoeffs[,c("InsLength", ResponseVarName)],
                                by = list(InsLBins, L1SingletonCoeffs$blnSelect), 
                                FUN = function(x) median(x, na.rm = T))
  }
  
  # Calculate standard error of response variable per LINE-1 length class
  StErrPerL <- aggregate(L1SingletonCoeffs[, ResponseVarName],
                             by = list(InsLBins), 
                             FUN = function(x) sqrt(var(x, na.rm = T) / length(x)))
  
  StErrPerLSelect  <- aggregate(L1SingletonCoeffs[,ResponseVarName],
                                    by = list(InsLBins, L1SingletonCoeffs$blnSelect), 
                                    FUN = function(x) sqrt(var(x, na.rm = T) / length(x)))
  
  # Difference between non-selected L1 and overall mean
  DiffNonSelectMean <- MeanPerLSelect[!MeanPerLSelect$Group.2, ResponseVarName] - 
    MeanVal
  
  # Difference between selected and non-selected L1 and overall mean
  SelectMatch  <- match(MeanPerLSelect$Group.1[MeanPerLSelect$Group.2],
                        MeanPerLSelect$Group.1[!MeanPerLSelect$Group.2])
  DiffSelectNonSelect <- MeanPerLSelect[MeanPerLSelect$Group.2, ResponseVarName] - 
    MeanPerLSelect[!MeanPerLSelect$Group.2, ResponseVarName][SelectMatch]
  
  # Test for correlation between difference between non-selected L1 and overall
  # mean and difference between selected and non-selected L1
  CorTestP <- sum(cor(DiffNonSelectMean[SelectMatch], DiffSelectNonSelect) >=
                    SampledCor) / NrSamples 
  
  # Put output in a list
  list(MeanVal = MeanVal, MeanPerL = MeanPerL, MeanPerLSelect = MeanPerLSelect,
       SampledVar = SampledVar, DiffNonSelectMean = DiffNonSelectMean[SelectMatch],
       DiffSelectNonSelect = DiffSelectNonSelect, CorTestP = CorTestP,
       StErrPerL = StErrPerL, StErrPerLSelect = StErrPerLSelect)
}  
  
SelPerInsList_Dist2GeneEnh <- SelectPatternPerInsL(ResponseVarName = "Dist2Enhancer")
SelPerInsList_Dist2GeneMean <- SelectPatternPerInsL(
  InsLBins = cut(L1SingletonCoeffs$InsLength, breaks = seq(0, 6300, 475)))
SelPerInsList_Dist2GeneMed  <- SelectPatternPerInsL(Method = "Median")
1 - SelPerInsList_Dist2GeneMean$CorTestP
SelPerInsList_Dist2GeneEnh$CorTestP
SelPerInsList_Dist2GeneMed$CorTestP
plot(SelPerInsList_Dist2GeneMed$DiffNonSelectMean,
     SelPerInsList_Dist2GeneMed$DiffSelectNonSelect)


# Plot mean distance to gene against insertion length
SelPerInsPlot <- function(SelPerInsList) {
  plot(SelPerInsList$MeanPerL$InsLength, 
       SelPerInsList$MeanPerL$Dist2Gene,
       ylim = c(0, 5*10^5), xlab = "L1 insertion length [bp]",
       ylab = "Mean distance to closest gene [Mb]",
       yaxt = "n")
  axis(2, at = 0:5*10^5, 0:5/10)
  QMat_Dist2Gene <- QMat(SelPerInsList$SampledVar, xVals  = 
                           SelPerInsList$MeanPerL$InsLength)
  points(SelPerInsList$MeanPerL$InsLength, 
         SelPerInsList$MeanPerL$Dist2Gene)
  points(SelPerInsList$MeanPerLSelect$InsLength[SelPerInsList$MeanPerLSelect$Group.2], 
         SelPerInsList$MeanPerLSelect$Dist2Gene[SelPerInsList$MeanPerLSelect$Group.2],
         pch = 16)
  lines(c(0, 10^5), rep(SelPerInsList$MeanVal, 2), lty = 2) 
  AddErrorBars(MidX = SelPerInsList$MeanPerL$InsLength, 
               MidY = SelPerInsList$MeanPerL$Dist2Gene,
               ErrorRange = SelPerInsList$StErrPerL$x, 
               TipWidth = 100)
  SelPerInsList$StErrPerL
  AddErrorBars(MidX = SelPerInsList$MeanPerLSelect$InsLength[
    SelPerInsList$MeanPerLSelect$Group.2], 
    MidY = SelPerInsList$MeanPerLSelect$Dist2Gene[
      SelPerInsList$MeanPerLSelect$Group.2],
    ErrorRange = SelPerInsList$StErrPerLSelect$x[
      SelPerInsList$StErrPerLSelect$Group.2], 
    TipWidth = 100)
  legend("topleft", bty = "n", legend = c("selected L1", "all L1"),
         pch = c(16, 1), cex = 0.75, y.intersp = 0.5)
  
}
SelPerInsPlot(SelPerInsList_Dist2GeneMean)
CreateDisplayPdf('D:/L1polymORF/Figures/Dist2GeneVsInsLength.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 5, width = 5)

##########################
#                        #
#    Save image          #
#                        #
##########################

# save.image(OutputPath)

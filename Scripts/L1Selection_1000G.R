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
L1FstPath       <- 'D:/L1polymORF/Data/L1all.weir.fst.weir.fst'
L1FstPath_pop   <- 'D:/L1polymORF/Data/L1all.weir.fst.Pops.weir.fst'

# Specify parameters
NrInfoCols <- 9

# Maximum L1 insertion length to be labeled a 'fragment'
MaxFragLength <- 5900

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

# Load previously generated objects
load(L1RefRangePath)
load(L1GRPath)

# Indicator for minimum frequency
blnAboveMinFreq <- L1_1000G_reduced$Frequency > 1/2504
log(10/2504)
# Read repeat table and subset to get only L1HS rows with fragment size below 
# MaxFragLength
RepeatTable      <- read.csv("D:/L1polymORF/Data/repeatsHg38_L1HS.csv")
RepeatTable_hg19 <- read.csv("D:/L1polymORF/Data/repeatsHg19_L1HS.csv")

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
L1RefGR <- GRanges(seqnames = RepeatTable$genoName,
                   ranges = IRanges(start = RepeatTable$genoStart,
                                    end = RepeatTable$genoEnd),
                   strand = RepeatTable$strand)
L1RefGRFull  <- L1RefGR[width(L1RefGR) > 6000]
L1FragmGR    <- L1RefGR[width(L1RefGR) < MaxFragLength]
L1RefGR_hg19 <- GRanges(seqnames = RepeatTable_hg19$genoName,
                        ranges = IRanges(start = RepeatTable_hg19$genoStart,
                                         end = RepeatTable_hg19$genoEnd),
                        strand = RepeatTable_hg19$strand)
L1RefGRFull_hg19 <- L1RefGR_hg19[width(L1RefGR_hg19) > 6000]
L1FragmGR_hg19   <- L1RefGR_hg19[width(L1RefGR_hg19) < MaxFragLength]

##########
# Process Fst Data
##########

# read Fst Data
L1Fst     <- read.delim(L1FstPath)
L1Fst_pop <- read.delim(L1FstPath_pop)

# Paste chromosome and position from Fst and 1000 Genome data
chrPos_Fst       <- paste(L1Fst$CHROM, L1Fst$POS, sep = "_")
chrPos_Fst_pop   <- paste(L1Fst_pop$CHROM, L1Fst_pop$POS, sep = "_")
chrPos_1000G     <- paste(L1_1000G$CHROM, L1_1000G$POS, sep = "_")

# Match Fst to 1000 genome data
ChrPosMatch       <- match(chrPos_1000G, chrPos_Fst)
ChrPosMatch_pop   <- match(chrPos_1000G, chrPos_Fst_pop)
L1Fst_matched     <- L1Fst[ChrPosMatch,]
L1Fst_matched_pop <- L1Fst_pop[ChrPosMatch_pop,]

##########
# Process L1 catalog
##########

cat("Matching 1000G L1 to L1 catalog\n")

# Read in table with Info about 1000 genome samples 
SampleInfo_1000Genome <- read.table(G1000SamplePath, as.is = T, header = T)

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

DistCat2_1000G <- Dist2Closest(L1CatalogGR, L1_1000G_GRList_hg38$LiftedRanges)

# Get indices of 1000 Genome and catalog elements that match
idx1000G <- nearest(L1CatalogGR, L1_1000G_GRList_hg38$LiftedRanges)
L1_1000G_GRList_hg38$idxUniqueMapped
L1CatalogMatch1000G <- L1CatalogL1Mapped[DistCat2_1000G < 100, ]
idx1000GMatchCat    <- L1_1000G_GRList_hg38$idxUniqueMapped[idx1000G[DistCat2_1000G < 100]]
L1CatalogGRMatched  <- L1CatalogGR[DistCat2_1000G < 100]
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

# List of files with domains and loops
DomainFiles <- list.files("D:/L1polymORF/Data/HiCData/", pattern = "domainlist.txt",
                          full.names = T)
DomainFiles <- DomainFiles[! DomainFiles %in% grep(".gz", DomainFiles, value = T)]
LoopFiles <- list.files("D:/L1polymORF/Data/HiCData/", pattern = "looplist.txt",
                        full.names = T)
LoopFiles <- LoopFiles[! LoopFiles %in% grep(".gz", LoopFiles, value = T)]

# Get list of domain and loop ranges
DomainGRList   <- lapply(DomainFiles, getLoopGRs)
LoopGRList     <- lapply(LoopFiles, getLoopGRs)
names(DomainGRList) <- sapply(DomainFiles, 
                              function(x) strsplit(x, "_")[[1]][2])
names(LoopGRList) <- sapply(LoopFiles, 
                            function(x) strsplit(x, "_")[[1]][2])

# Summarize domain ranges (union or intersect)
DomainGR_Intersect <- DomainGRList[[1]]$AllLoops
LoopGR_Intersect   <- LoopGRList[[1]]$AllLoops
for (i in 2:length(LoopGRList)){
  DomainGR_Intersect <- intersect(DomainGR_Intersect, DomainGRList[[i]]$AllLoops)
  LoopGR_Intersect   <- intersect(LoopGR_Intersect, LoopGRList[[i]]$AllLoops)
}


##########################################
#                                        #
#           Fst                          #
#                                        #
##########################################

# Fst as function of insertion length and 
LMF <- lm(L1Fst_matched$WEIR_AND_COCKERHAM_FST ~ L1_1000G$InsLength + blnFull)
summary(LMF)
LMF_minFreq <- lm(L1Fst_matched$WEIR_AND_COCKERHAM_FST[blnAboveMinFreq] ~ 
                    L1_1000G$InsLength[blnAboveMinFreq] + blnFull[blnAboveMinFreq])
summary(LMF_minFreq)

# Summarize Fst per insertion length 
InsLengthCut <- cut(L1_1000G_reduced$InsLength, breaks = seq(0, 6500, 500))
FstPerInsL   <- aggregate(cbind(L1Fst_matched$WEIR_AND_COCKERHAM_FST,
                                     L1Fst_matched_pop$WEIR_AND_COCKERHAM_FST,
                                     L1_1000G_reduced$InsLength), 
                               by = list(InsLengthCut), 
                               FUN = function(x) mean(x, na.rm = T))
FstPerInsL_Var  <- aggregate(cbind(L1Fst_matched$WEIR_AND_COCKERHAM_FST,
                                       L1Fst_matched_pop$WEIR_AND_COCKERHAM_FST,
                                       L1_1000G_reduced$InsLength), 
                                 by = list(InsLengthCut), 
                                 FUN = function(x) var(x, na.rm = T))

# Summarize Fst per insertion length for minimum frequencies
InsLengthCut <- cut(L1_1000G_reduced$InsLength, breaks = seq(0, 6500, 500))
FstPerInsL_MinF   <- aggregate(cbind(L1Fst_matched$WEIR_AND_COCKERHAM_FST[blnAboveMinFreq],
                              L1Fst_matched_pop$WEIR_AND_COCKERHAM_FST[blnAboveMinFreq],
                              L1_1000G_reduced$InsLength[blnAboveMinFreq]), 
                        by = list(InsLengthCut[blnAboveMinFreq]), 
          FUN = function(x) mean(x, na.rm = T))
FstPerInsL_MinFVar  <- aggregate(cbind(L1Fst_matched$WEIR_AND_COCKERHAM_FST[blnAboveMinFreq],
                              L1Fst_matched_pop$WEIR_AND_COCKERHAM_FST[blnAboveMinFreq],
                              L1_1000G_reduced$InsLength[blnAboveMinFreq]), 
                        by = list(InsLengthCut[blnAboveMinFreq]), 
                        FUN = function(x) var(x, na.rm = T))

plot(FstPerInsL_MinF$V3, FstPerInsL_MinF$V1, ylim = c(0, 0.15), xlab = "Insertion length",
     ylab = "Fst")
points(FstPerInsL$V3, FstPerInsL$V1, col = "red")
AddErrorBars(MidX = FstPerInsL_MinF$V3, MidY = FstPerInsL_MinF$V1,
             ErrorRange = sqrt(FstPerInsL_MinFVar$V1 / 
                                 table(InsLengthCut[blnAboveMinFreq])),
             TipWidth = 100)
plot(FstPerInsL$V3, FstPerInsL$V2)

##########################################
#                                        #
#     Compare Fst and 
#     distance to genes among fragments  #
#                                        #
##########################################

# Get ranges of genes
GRgenes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Distances to closest gene
Dist2ClosestGene <- Dist2Closest(L1_1000G_GR_hg19, GRgenes_hg19) 
blnFull <- L1_1000G_reduced$InsLength >= 6000

# Correlation between distance and frequency
LMF <- lm(L1_1000G_reduced$Frequency[blnFull] ~ Dist2ClosestGene[blnFull])
summary(LMF)
cor.test(L1_1000G_reduced$Frequency[blnFull], Dist2ClosestGene[blnFull], method = "kendall")


cor.test(L1Fst_matched$WEIR_AND_COCKERHAM_FST[blnFull], 
         Dist2ClosestGene[blnFull], method = "kendall")
plot(L1Fst_matched$WEIR_AND_COCKERHAM_FST[blnFull] ~ Dist2ClosestGene[blnFull])

DistCut <- cut(Dist2ClosestGene, breaks = seq(0, 1.5*10^6, 10^5))
FstPerDist <- aggregate(cbind(L1Fst_matched$WEIR_AND_COCKERHAM_FST[blnFull],
                              L1Fst_matched_pop$WEIR_AND_COCKERHAM_FST[blnFull],
                              Dist2ClosestGene[blnFull]), by = list(DistCut[blnFull]), 
                        FUN = function(x) mean(x, na.rm = T))
plot(FstPerDist$V3, FstPerDist$V1)


# Distances to closest loop
Dist2ClosestLoop <- Dist2Closest(L1_1000G_GR_hg19, LoopGR_Intersect) 
blnDistLoopBelow <- Dist2ClosestLoop < 7*10^4

LMF <- lm(L1Fst_matched$WEIR_AND_COCKERHAM_FST[blnFull] ~ 
            Dist2ClosestLoop[blnFull] + L1_1000G_reduced$Frequency[blnFull] +
            Dist2ClosestLoop[blnFull]:L1_1000G_reduced$Frequency[blnFull])
LMF <- lm(L1Fst_matched$WEIR_AND_COCKERHAM_FST[blnAboveMinFreq & blnFull] ~ 
            Dist2ClosestLoop[blnAboveMinFreq & blnFull])
LMF2 <- lm(L1Fst_matched$WEIR_AND_COCKERHAM_FST[!blnFull] ~ Dist2ClosestLoop[!blnFull])
LMF3 <- lm(L1Fst_matched$WEIR_AND_COCKERHAM_FST[blnAboveMinFreq] ~ 
             Dist2ClosestLoop[blnAboveMinFreq] + blnFull[blnAboveMinFreq]
           + Dist2ClosestLoop[blnAboveMinFreq]:blnFull[blnAboveMinFreq])
summary(LMF)
summary(LMF2)
summary(LMF3)
plot(L1Fst_matched$WEIR_AND_COCKERHAM_FST[blnFull] ~ Dist2ClosestLoop[blnFull],
     xlim = c(0, 10^6))
plot(L1Fst_matched$WEIR_AND_COCKERHAM_FST[blnAboveMinFreq & blnFull] ~ 
          Dist2ClosestLoop[blnAboveMinFreq & blnFull])
plot(L1Fst_matched$WEIR_AND_COCKERHAM_FST[!blnFull] ~ Dist2ClosestLoop[!blnFull],
     xlim = c(0, 10^6))

# Distances to closest domain
Dist2ClosestDomain <- Dist2Closest(L1_1000G_GR_hg19, DomainGR_Intersect) 
LMF <- lm(L1Fst_matched$WEIR_AND_COCKERHAM_FST[idx1000GMatchCat] ~ Dist2ClosestDomain[idx1000GMatchCat])
summary(LMF)

# Fst vs activity
LMF <- lm(L1Fst_matched$WEIR_AND_COCKERHAM_FST[idx1000GMatchCat] ~ 
            L1CatalogMatch1000G$ActivityNum)
summary(LMF)
LMF <- lm(L1Fst_matched$WEIR_AND_COCKERHAM_FST[idx1000GMatchCat] ~ 
            L1CatalogMatch1000G$ActivityNum + Dist2ClosestDomain[idx1000GMatchCat])
summary(LMF)

plot(log(L1_1000G_reduced$Frequency), L1Fst_matched$WEIR_AND_COCKERHAM_FST)
plot(L1_1000G_reduced$Frequency, L1Fst_matched$WEIR_AND_COCKERHAM_FST)

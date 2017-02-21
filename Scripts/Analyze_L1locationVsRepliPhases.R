# The following script reads in a a bed file of replication phases and determines 
# whether full-length L1 are overrepresented in a particular phase

library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)
library(Hotelling)

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.r')

############
#  Process L1 data
############

# Read repeat table with L1HS in hg19
RepeatTable <- read.table("D:/L1polymORF/Data/repeatsHg19_L1HS", header = T)

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
L1RefGR <- GRanges(seqnames = RepeatTable$genoName,
                   ranges = IRanges(start = RepeatTable$genoStart,
                                    end = RepeatTable$genoEnd),
                   strand = RepeatTable$strand)

# Path to L1 catalogue file 
L1Catalogue <- read.csv("D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv", as.is = T)


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
L1CatGR_hg19 <- liftOver(L1CatalogGR, 
                         chain = import.chain("D:/L1polymORF/Data/hg38ToHg19.over.chain"))
NrRanges <- sapply(L1CatGR_hg19, length)
idxUniqueMapped <- NrRanges == 1
L1CatGR_hg19 <- unlist(L1CatGR_hg19[idxUniqueMapped])
L1CatGR_hg19 <- L1CatGR_hg19[width(L1CatGR_hg19) > 6000]

# Find overlap with fragments and catalog elements
L1RefGR_fragm    <- L1RefGR[width(L1RefGR) < 6000]

# Get combined genomic ranges and create a bed file
L1GR_combined <- c(L1CatGR_hg19, L1RefGR_fragm)
L1GR_combinedDF <- as.data.frame(L1GR_combined)
write.table(L1GR_combinedDF[1:1000,1:3], "D:/L1polymORF/Data/L1GRtable1",
            quote = F, row.names = F, col.names = F)
write.table(L1GR_combinedDF[1001:nrow(L1GR_combinedDF),1:3], "D:/L1polymORF/Data/L1GRtable2",
            quote = F, row.names = F, col.names = F)
write.table(L1GR_combinedDF[,1:3], "D:/L1polymORF/Data/L1GRtableFull",
            quote = F, row.names = F, col.names = F)

############
#  Process replication data
############

# Read table with open compartment data
RepliTableFiles <- list.files("D:/L1polymORF/Data/UWrepliSeqData", 
   pattern = "UWrepliSeqGM12878_L1GRall", full.names = T)
# RepliTableFiles <- list.files("D:/L1polymORF/Data/", 
#                               pattern = "UWrepliSeqGM12878_L1GR", full.names = T)
# RepliTableFiles <- RepliTableFiles[-1]

# Function to make a replication phase table
makeRepliTable <- function(FileName) {
  NameSplt   <- strsplit(FileName, "_")[[1]]
  RepliPhase <- NameSplt[length(NameSplt)]
  RepliLines <- readLines(FileName)
  idxNewL1 <- grep("variableStep", RepliLines)
  ID      <- 1:length(idxNewL1)
  NrRepl   <- c(idxNewL1[-1], length(RepliLines) + 1) - idxNewL1 - 1
  Chroms   <- sapply(RepliLines[idxNewL1], function(x) {
    SpltLine <- strsplit(x, " ")[[1]]
    strsplit(SpltLine[2], "=")[[1]][2]
  })
  startV               <- rep(NA, length(idxNewL1))
  startLines           <- RepliLines[idxNewL1[NrRepl > 0] + 1]
  StartTable           <- read.delim(text = startLines, header = F)
  endV                 <- rep(NA, length(idxNewL1))
  endLines             <- RepliLines[idxNewL1[NrRepl > 0] + NrRepl[NrRepl > 0]]
  endTable             <- read.delim(text = endLines, header = F)
  startV[NrRepl > 0]   <- StartTable[,1]
  endV[NrRepl > 0]     <- endTable[,1]
  RepliTable           <- read.delim(text = RepliLines[-c(1, idxNewL1)], header = F)
  colnames(RepliTable) <- c("startLocal", "PhasePercent")
  RepliTable$chrom     <- rep(Chroms, NrRepl)
  RepliTable$start     <- rep(startV, NrRepl)
  RepliTable$end       <- rep(endV, NrRepl)
  RepliTable$ID        <- rep(ID, NrRepl)
  RepliAgg             <- aggregate(cbind(start, end, PhasePercent) ~ ID, 
                                    data = RepliTable, FUN = mean)
  RepliAgg$chrom <- Chroms[NrRepl > 0]
  colnames(RepliAgg)[colnames(RepliAgg) == "PhasePercent"] <- RepliPhase
  return(RepliAgg[ ,colnames(RepliAgg) != "ID"])
}

# Create an aggregated replication table
RepliAgg <- makeRepliTable(RepliTableFiles[1])
for (FileName in RepliTableFiles[-1]){
  NewTable <- makeRepliTable(FileName)
  RepliAgg <- merge(RepliAgg, NewTable)
}

# Turn open chromatin table into genomic ranges
RepliGR          <- makeGRangesFromDataFrame(RepliAgg)
idxOverlap_fragm <- findOverlaps(RepliGR, L1RefGR_fragm)
idxOverlap_Cat   <- findOverlaps(RepliGR, L1CatGR_hg19)
idxOverlap_fragm@from
max(idxOverlap_Cat@to)

# Get 2 replication matrices, 1 for fragments and one for full L1
OrderedPhases <- c("G1b", "S1", "S2", "S3", "S4", "G2")
#OrderedPhases <- c("S1", "S2", "S3", "S4", "G2")
RepliFragm <- RepliAgg[idxOverlap_fragm@from, OrderedPhases]
RepliFull  <- RepliAgg[idxOverlap_Cat@from, OrderedPhases]
StdErr_Fragm <- apply(RepliFragm, 2, FUN = function(x) sqrt(var(x) / length(x)))
StdErr_Full  <- apply(RepliFull, 2, FUN = function(x) sqrt(var(x) / length(x)))
plot(colMeans(RepliFragm), type = "l", col = "red", xaxt = "n", 
     ylim = c(5, 27), xlab = "Replication phase", ylab = "Percentage")
axis(1, at = seq_along(OrderedPhases), labels = OrderedPhases)
lines(colMeans(RepliFull), col = "blue")
AddErrorBars(MidX = seq_along(OrderedPhases), MidY = colMeans(RepliFragm),
             ErrorRange = StdErr_Fragm, TipWidth = 0.02, Col = "red")
AddErrorBars(MidX = seq_along(OrderedPhases), MidY = colMeans(RepliFull),
             ErrorRange = StdErr_Full, TipWidth = 0.02, Col = "blue")
FitDataRepliFragm <- matrix(nrow = nrow(RepliFragm) * length(OrderedPhases), ncol = 2)
FitDataRepliFull  <- matrix(nrow = nrow(RepliFull) * length(OrderedPhases), ncol = 2)
for (i in 1:length(OrderedPhases)){
  idxFragm <- ((i - 1) * nrow(RepliFragm) + 1):(i * nrow(RepliFragm))
  idxFull  <- ((i - 1) * nrow(RepliFull) + 1):(i * nrow(RepliFull))
  FitDataRepliFragm[idxFragm, 1] <- rep(i, nrow(RepliFragm))
  FitDataRepliFull[idxFull, 1]  <- rep(i, nrow(RepliFull))
  FitDataRepliFragm[idxFragm, 2] <- RepliFragm[,i]
  FitDataRepliFull[idxFull, 2]  <- RepliFull[,i]
}
model_Fragm <- lm(FitDataRepliFragm[,2] ~ poly(FitDataRepliFragm[,1], 2))
model_Full  <- lm(FitDataRepliFull[,2] ~ poly(FitDataRepliFull[,1], 2))


# Group chromatin data by overlap
G2        <- RepliAgg$G2
S2        <- RepliAgg$S2
WeightedPhase <- as.matrix(RepliAgg[, OrderedPhases]) %*% 1:6
RepliByOverlap <- data.frame(
  Percent_S2 = c(G2[idxOverlap_fragm@from], G2[idxOverlap_Cat@from]),
  Percent_G2 = c(G2[idxOverlap_fragm@from], G2[idxOverlap_Cat@from]),
  WgtPhase   = c(WeightedPhase[idxOverlap_fragm@from], 
                 WeightedPhase[idxOverlap_Cat@from]),
  L1Type     = c(rep("Fragm", length(idxOverlap_fragm@from)),
                 rep("Full", length(idxOverlap_Cat@from))))
boxplot(RepliByOverlap$Percent_G2 ~ RepliByOverlap$L1Type)
boxplot(RepliByOverlap$WgtPhase ~ RepliByOverlap$L1Type)
t.test(Percent_S2 ~ L1Type, data = RepliByOverlap)
t.test(Percent_G2 ~ L1Type, data = RepliByOverlap)
t.test(WgtPhase ~ L1Type, data = RepliByOverlap)
hist(RepliByOverlap$WgtPhase[RepliByOverlap$L1Type == "Fragm"])
hist(RepliByOverlap$WgtPhase[RepliByOverlap$L1Type == "Full"])


# Distinguish full-length and fragment
idxOverlap <- findOverlaps(RepliGR, L1RefGR)
blnFullL1  <- width(L1RefGR[idxOverlap@to]) > 6000
boxplot(WeightedPhase[idxOverlap@from] ~ blnFullL1)
t.test(WeightedPhase[idxOverlap@from] ~ blnFullL1)

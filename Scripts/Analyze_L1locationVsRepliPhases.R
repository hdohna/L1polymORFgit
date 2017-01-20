# The following script reads in a a bed file of replication phases and determines 
# whether full-length L1 are overrepresented in a particular phase

library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)
library(Hotelling)

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

# Find overlap with fragments and catalog elements
L1RefGR_fragm    <- L1RefGR[width(L1RefGR) < 5000]

# Get combined genomic ranges and create a bed file
L1GR_combined <- c(L1CatGR_hg19, L1RefGR_fragm)
L1GR_combinedDF <- as.data.frame(L1GR_combined)
write.table(L1GR_combinedDF[1:1000,1:3], "D:/L1polymORF/Data/L1GRtable",
            quote = F, row.names = F, col.names = F)

############
#  Process replication data
############

# Read table with open compartment data
RepliTableFiles <- list.files("D:/L1polymORF/Data/", 
   pattern = "UWrepliSeqGM12878_L1GR", full.names = T)
RepliTable <- read.table(RepliTableFiles[1], skip = 1)

FileName   <- RepliTableFiles[2]
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
RepliAgg <- makeRepliTable(RepliTableFiles[2])
for (FileName in RepliTableFiles[-c(1:2)]){
  NewTable <- makeRepliTable(FileName)
  RepliAgg <- merge(RepliAgg, NewTable)
}

# Turn open chromatin table into genomic ranges
RepliGR          <- makeGRangesFromDataFrame(RepliAgg)
idxOverlap_fragm <- findOverlaps(RepliGR, L1RefGR_fragm)
idxOverlap_Cat   <- findOverlaps(RepliGR, L1CatGR_hg19)

# Group chromatin data by overlap
G2        <- RepliAgg$G2
WeightedPhase <- as.matrix(RepliAgg[, c("S1", "S2", "S3", "S4", "G2")]) %*% c(1:5)
RepliByOverlap <- data.frame(
  Percent_G2 = c(G2[idxOverlap_fragm@queryHits], G2[idxOverlap_Cat@queryHits]),
  WgtPhase   = c(WeightedPhase[idxOverlap_fragm@queryHits], 
                 WeightedPhase[idxOverlap_Cat@queryHits]),
  L1Type     = c(rep("Fragm", length(idxOverlap_fragm@queryHits)),
                 rep("Full", length(idxOverlap_Cat@queryHits))))
boxplot(RepliByOverlap$Percent_G2 ~ RepliByOverlap$L1Type)
boxplot(RepliByOverlap$WgtPhase ~ RepliByOverlap$L1Type)
t.test(Percent_G2 ~ L1Type, data = RepliByOverlap)
t.test(WgtPhase ~ L1Type, data = RepliByOverlap)

HT <- hotelling.test(RepliAgg[idxOverlap_fragm@queryHits, c("G2", "S1", "S2", "S3", "S4")],
               RepliAgg[idxOverlap_Cat@queryHits, c("G2", "S1", "S2", "S3", "S4")])
HT

# Distinguish full-length and fragment
idxOverlap <- findOverlaps(RepliGR, L1RefGR)
blnFullL1  <- width(L1RefGR[idxOverlap@subjectHits]) > 6000
boxplot(WeightedPhase[idxOverlap@queryHits] ~ blnFullL1)
t.test(WeightedPhase[idxOverlap@queryHits] ~ blnFullL1)

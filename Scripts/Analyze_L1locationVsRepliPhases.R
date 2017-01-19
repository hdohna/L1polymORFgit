# The following script reads in a a bed file of replication phases and determines 
# whether full-length L1 are overrepresented in a particular phase

library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)

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


############
#  Process replication data
############

# Read table with open compartment data
RepliTableFiles <- list.files("D:/L1polymORF/Data/", 
   pattern = "UWrepliSeqPerL1", full.names = T)
RepliTable <- read.table(RepliTableFiles[3], skip = 1)
for (FileName in RepliTableFiles[-1]){
  NewTable <- read.table(FileName, skip = 1)
  RepliTable <- merge(RepliTable, NewTable)
}
RepliTable <- import.bed("D:/L1polymORF/Data/UWrepliSeq.bed")
RepliTable <- readLines("D:/L1polymORF/Data/UWrepliSeq.bed")
RepliTable[nchar(RepliTable) > 20]
RepliTable[1:10]
unique(SubCTable$chrom)

# Turn open chromatin table into genomic ranges
SubCGR <- makeGRangesFromDataFrame(SubCTable)

# Find overlap with fragments and catalog elements
L1RefGR_fragm    <- L1RefGR[width(L1RefGR) < 5000]
idxOverlap_fragm <- findOverlaps(SubCGR, L1RefGR_fragm)
idxOverlap_Cat   <- findOverlaps(SubCGR, L1CatGR_hg19)

# Get combined genomic ranges and create a bed file
L1GR_combined <- c(L1CatGR_hg19, L1RefGR_fragm)
L1GR_combinedDF <- as.data.frame(L1GR_combined)
write.table(L1GR_combinedDF[1:1000,1:3], "D:/L1polymORF/Data/L1GRtable",
            quote = F, row.names = F, col.names = F)

# Group chromatin data by overlap
SubC_char <- as.character(SubCTable$SubC)
AB <- substr(SubC_char, 1, 1)
SubCByOverlap <- data.frame(
  AB = c(AB[idxOverlap_fragm@queryHits], AB[idxOverlap_Cat@queryHits]),
  SubC = c(SubC_char[idxOverlap_fragm@queryHits],
                                     SubC_char[idxOverlap_Cat@queryHits]),
  L1Type = c(rep("Fragm", length(idxOverlap_fragm@queryHits)),
                                      rep("Full", length(idxOverlap_Cat@queryHits))))
ConTable <- table(SubCByOverlap$SubC, SubCByOverlap$L1Type)
ConTable[,2] / rowSums(ConTable)
chisq.test(SubCByOverlap$SubC, SubCByOverlap$L1Type)
chisq.test(SubCByOverlap$AB, SubCByOverlap$L1Type)
fisher.test(SubCByOverlap$SubC, SubCByOverlap$L1Type)
fisher.test(SubCByOverlap$AB, SubCByOverlap$L1Type)

# Distinguish full-length and fragment
idxOverlap <- findOverlaps(HiCGR, L1RefGR_hg18)
idxOverlap@subjectHits

blnFullL1 <- width(L1RefGR_hg18[idxOverlap@subjectHits]) > 6000
boxplot(HiCData$Eij[idxOverlap@queryHits] ~ blnFullL1)
aggPerWindow <- aggregate(blnFullL1 ~ idxOverlap@queryHits, FUN = mean)
aggPerWindow$Eij <- HiCData$Eij[aggPerWindow$`idxOverlap@queryHits`]
plot(aggPerWindow$Eij, aggPerWindow$blnFullL1)

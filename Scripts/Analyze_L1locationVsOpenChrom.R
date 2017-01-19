# The following script reads in a Hi-C matrix and determines well-connected stretchs

library(GenomicRanges)
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
                         chain = import.chain(ChainFile38To19))
NrRanges <- sapply(L1CatGR_hg19, length)
idxUniqueMapped <- NrRanges == 1
L1CatGR_hg19 <- unlist(L1CatGR_hg19[idxUniqueMapped])


############
#  Process chromatin data
############

# Read table with open chromatin data
ChromTable <- read.table("D:/L1polymORF/Data/openChromSynthGm12878", header = T)

# Turn open chromatin table into genomic ranges
ChromGR <- GRanges(seqnames = ChromTable$chrom,
                   ranges = IRanges(start = ChromTable$chromStart,
                                    end = ChromTable$chromEnd))

# Find overlap with fragments and catalog elements
L1RefGR_fragm    <- L1RefGR[width(L1RefGR) < 5000]
idxOverlap_fragm <- findOverlaps(ChromGR, L1RefGR_fragm)
idxOverlap_Cat   <- findOverlaps(ChromGR, L1CatGR_hg19)

# Group chromatin data by overlap
ChromByOverlap <- data.frame(HiC = c(HiCData$Eij[idxOverlap_fragm@queryHits],
                                   HiCData$Eij[idxOverlap_Cat@queryHits]),
                           L1Type = c(rep("Fragm", length(idxOverlap_fragm@queryHits)),
                                      rep("Full", length(idxOverlap_Cat@queryHits))))
boxplot(HiCByOverlap$HiC ~ HiCByOverlap$L1Type)

# Distinguish full-length and fragment
idxOverlap <- findOverlaps(HiCGR, L1RefGR_hg18)
idxOverlap@subjectHits

blnFullL1 <- width(L1RefGR_hg18[idxOverlap@subjectHits]) > 6000
boxplot(HiCData$Eij[idxOverlap@queryHits] ~ blnFullL1)
aggPerWindow <- aggregate(blnFullL1 ~ idxOverlap@queryHits, FUN = mean)
aggPerWindow$Eij <- HiCData$Eij[aggPerWindow$`idxOverlap@queryHits`]
plot(aggPerWindow$Eij, aggPerWindow$blnFullL1)

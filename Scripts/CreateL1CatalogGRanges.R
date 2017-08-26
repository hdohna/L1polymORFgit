# The script below reads in L1 catalog, creates genomic ranges and saves them

# Load packages
library(rtracklayer)

# Path to L1 catalogue file (Created in script AddColumns2L1Catalog.R)
L1CataloguePath   <- "D:/L1polymORF/Data/L1CatalogExtended.csv"
L1CatalogGROutput <- "D:/L1polymORF/Data/L1CatalogGRanges.RData"

# Read in table with known L1 
L1Catalogue <- read.csv(L1CataloguePath, as.is = T)
L1Catalogue$Allele[is.na(L1Catalogue$Allele)] <- 1

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
blnInRef1         <- (L1Catalogue$end_HG38 - L1Catalogue$start_HG38) > 6000 
L1CatalogL1Mapped <- L1Catalogue[blnL1Mapped & blnAllele1,]

# Lift catalog ranges to hg19
L1LiftoverList <- LiftoverL1Catalog(L1CatalogL1Mapped,  
                                    ChainFilePath = "D:/L1polymORF/Data/hg38ToHg19.over.chain")
L1CatalogGR_hg19 <- L1LiftoverList$GRCatalogue_hg19

# Create genomic ranges for catalog L1
L1CatalogGR <- GRanges(seqnames = L1CatalogL1Mapped$Chromosome,
                       ranges = IRanges(start = pmin(L1CatalogL1Mapped$start_HG38,
                                                     L1CatalogL1Mapped$end_HG38),
                                        end = pmax(L1CatalogL1Mapped$start_HG38,
                                                   L1CatalogL1Mapped$end_HG38)),
                       strand = L1CatalogL1Mapped$strand_L1toRef)

save.image(L1CatalogGROutput)

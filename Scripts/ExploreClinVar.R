# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

library(GenomicRanges)

##########
# Process ClinVar data
##########

# Read vcf file
ClinVarVcf <- ReadVCF("D:/L1polymORF/Data/clinvar_20170404_hg38.vcf")

# Create genomic ranges of clinvar
ClinVarVcf$chromosome <- paste("chr", ClinVarVcf$X.CHROM, sep = "")
ClinVarGR <- makeGRangesFromDataFrame(ClinVarVcf, start.field="POS",
                                      end.field=c("POS"))

##########
# Process L1 catalog
##########

# Read in table with known L1 
L1Catalogue <- read.csv("D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_16-33-31_2016.csv", as.is = T)
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

# Check whether any catalog L1 overlaps with Clinvar
ClinVarWideGR <- resize(ClinVarGR, 10000, fix = "center")
Overlaps      <- findOverlaps(L1CatalogGR, ClinVarWideGR)
L1CatalogClin <- L1CatalogL1Mapped[unique(Overlaps@from), ]
ClinVarVcf[unique(Overlaps@to), ]
sum(overlapsAny(L1CatalogGR, ClinVarWideGR))
sum(overlapsAny(ClinVarWideGR, L1CatalogGR))

# Shift L1s randomly and count how many clinVars they are associated with
OverlapSample <- sapply(1:100, function(x){
  ShiftVals <- sample((-500000):500000, length(L1CatalogGR))
  L1CatalogGR_shifted <- shift(L1CatalogGR, ShiftVals)
  sum(overlapsAny(L1CatalogGR_shifted, ClinVarWideGR))
})
mean(OverlapSample)

OverlapSampleClin <- sapply(1:100, function(x){
  ShiftVals <- sample((-500000):500000, length(L1CatalogGR))
  L1CatalogGR_shifted <- shift(L1CatalogGR, ShiftVals)
  sum(overlapsAny(ClinVarWideGR, L1CatalogGR_shifted))
})
mean(OverlapSampleClin)
hist(OverlapSampleClin, breaks = seq(0, 800, 20))

mean(OverlapSampleClin <= sum(overlapsAny(ClinVarWideGR, L1CatalogGR)))


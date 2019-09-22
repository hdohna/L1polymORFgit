# Load necessary packages
library(rtracklayer)
library(GenomicRanges)

# Read in tables by Beck et al. (2010)
BeckTable1 <- read.table("D:/OneDrive - American University of Beirut/L1polymORF/Data/Beck2020Table1.txt",
           header = T, as.is = T)
BeckTableS5 <- read.table("D:/OneDrive - American University of Beirut/L1polymORF/Data/Beck2010TableS5",
           header = T, as.is = T)
BeckTableS5$FosmidNr <- substr(BeckTableS5$Fosmid, 1, 1)
BeckTable <- merge(BeckTable1, BeckTableS5, all = T)

# The following script analyzes basic quantities of the MELT2.0 data by Gardener et al. (2017)
MEInsCall <- read.table("D:/OneDrive - American University of Beirut/L1polymORF/Data/nstd144.GRCh37.variant_call.vcf",
                        as.is = T,
                        col.names = c("Chrom", "Pos", "ID", "Alt", "Type", "V6", 
                                      "V7", "Info"))
MEInsCall <- MEInsCall[MEInsCall$Type == "<INS:ME:LINE1>",]
grep("L1Ta1", MEInsCall$Info)
MEInsCall$Info[1]

# Extract allele frequency from info column
GetAF <- function(x){
  xSplit <- strsplit(x, ";")[[1]]
  AFch   <- strsplit(xSplit[length(xSplit)], "=")[[1]][2]
  as.numeric(AFch)
}
GetLength <- function(x){
  xSplit <- strsplit(x, ";")[[1]]
  LengthCh   <- strsplit(xSplit[grep("SVLEN=", xSplit)], "=")[[1]][2]
  as.numeric(LengthCh)
}

# Add columns necessary for analysis 
MEInsCall$AF <- sapply(MEInsCall$Info, GetAF)
MEInsCall <- MEInsCall[!is.na(MEInsCall$AF), ]
MEInsCall$L1width <- sapply(MEInsCall$Info, GetLength)
MEInsCall$SampleSize <- 1/min(MEInsCall$AF) 
sum(MEInsCall$L1width >= 6000, na.rm = T)

# MEInsCall$SampleSize <- 2 * MEInsSamplesize
MEInsCall$Freq <- ceiling(MEInsCall$SampleSize * MEInsCall$AF) # TODO: Figure out why not integers!
MEInsCall$blnFull <- MEInsCall$L1width >= 6000

# Create genomic ranges for MEInsCall
MEInsCall$ChromName <- paste("chr", MEInsCall$Chrom, sep = "") 
GR_MEInsCall <- makeGRangesFromDataFrame(MEInsCall, seqnames.field = "ChromName",
                                         start.field = "Pos",
                                         end.field = "Pos")

# Read in vcf file with MELT deletion calls
MEDelCall <- ReadVCF("D:/OneDrive - American University of Beirut/L1polymORF/Data/DEL.final_comp.vcf")
MEDelCall$chromosome <- paste("chr", MEDelCall$X.CHROM, sep = "")
MEDel_GR  <- makeGRangesFromDataFrame(df = MEDelCall,
                                      start.field = "POS",
                                      end.field = "POS")
colnames(MEDelCall)
width(MEDel_GR)

# function to get numeric genotype
GetNumericGenotype <- function(x){
  Split1 <- strsplit(x, ":")[[1]][1]
  Split2 <- strsplit(Split1, "/")[[1]]
  sum(as.numeric(Split2))
}

# Get numeric genotype of all reference L1 deletions
GTCols <- grep("L1Filtered", colnames(MEDelCall))
L1RefNumGen <- 2 - sapply(GTCols, function(x){
  sapply(1:nrow(MEDelCall), function(y) GetNumericGenotype(MEDelCall[y,x]))
})

# Add columns for frequency and sample size
MEDelCall$Freq       <- rowSums(L1RefNumGen, na.rm = T)
MEDelCall$SampleSize <- apply(L1RefNumGen, 1, function(x) 2*sum(!is.na(x)))
MEDelCall$L1width <- sapply(MEDelCall$INFO, GetLength)
sum(MEDelCall$L1width >= 6000)
# Read in L1 catalog
L1Catalog <- read.csv("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1CatalogExtended.csv", as.is = T)
blnCoord <- (!is.na(L1Catalog$start_HG38)) & (!is.na(L1Catalog$end_HG38))
L1Catalog <- L1Catalog[blnCoord, ]
sum(L1Catalog$blnInRef & L1Catalog$ActivityNum > 0, na.rm = T)
sum(L1Catalog$Allele == 1)
sum(L1Catalog$Reference %in% c("Brouha2003", "Beck2010"))
sum(L1Catalog$Reference %in% c("Brouha2003", "Beck2010") & L1Catalog$ActivityNum > 0,
    na.rm = T)
LiftOverList <- LiftoverL1Catalog(L1Catalog,
                    ChainFilePath = "D:/OneDrive - American University of Beirut/L1polymORF/Data/hg38ToHg19.over.chain")
GRCatalogue_hg19 <- LiftOverList$GRCatalogue_hg19
GRCatalogue_hg19_large <- resize(GRCatalogue_hg19, 200, fix = "center")
L1Catalog_hg19   <- LiftOverList$L1CatalogWithHG19

# Check what proportion of non-reference and Beck L1s overlap with 1000G L1s
GRL1Beck_hg19 <- GRCatalogue_hg19[L1Catalog_hg19$Coriell_ID %in% SampleColumns]
GRL1Beck_hg19_large <- resize(GRL1Beck_hg19, 1000, fix = "center")

mean(overlapsAny(GRL1Beck_hg19_large, L1_1000G_GR_hg19))
sum(overlapsAny(GRL1Beck_hg19_large, L1_1000G_GR_hg19))
mean(overlapsAny(GRL1Beck_hg19_large, GR_MEInsCall))

MEInsCall$L1width[overlapsAny(GR_MEInsCall, GRL1Beck_hg19_large)]
L1_1000G$InsLength[overlapsAny(L1_1000G_GR_hg19, GRL1Beck_hg19_large)]

# Check for each coriell ID whether L1s were detected in 1000 Genome
CorrielIDs <- unique(L1Catalog_hg19$Coriell_ID)
CorrielIDs <- CorrielIDs[!is.na(CorrielIDs)]
L1Detection <- sapply(CorrielIDs[CorrielIDs %in% SampleColumns], function(x){
  GRCat   <- GRCatalogue_hg19_large[which(L1Catalog_hg19$Coriell_ID == x)]
  GR1000G <- L1_1000G_GR_hg19[which(L1_1000G[,x] > 0)]
  GR1000GFull <- L1_1000G_GR_hg19[which(L1_1000G[,x] > 0 & L1_1000G$InsLength >= 5900)]
  Nr1000G <- sum(L1_1000G[,x] > 0 & L1_1000G$InsLength >= 6000, na.rm = T)
  c(NrDetected = sum(overlapsAny(GRCat, GR1000G)), 
    NrPresentBeck = length(GRCat),
    NrPresent1KG  = Nr1000G,
    PropDetectedBeck = mean(overlapsAny(GRCat, GR1000G)),
  PropDetected1KG = mean(overlapsAny(GR1000GFull, GRCat)))
})

PropDetected <- sum(L1Detection["NrDetected",]) / sum(L1Detection["NrPresentBeck",]) 
sum(MEInsCall$L1width >= 6000, na.rm = T) / mean(overlapsAny(GRL1Beck_hg19_large, GR_MEInsCall))

idxFull <- which(L1_1000G$InsLength >= 6000)
mean(sapply(SampleColumns, function(x) sum(L1_1000G[idxFull,x] > 0))) / mean(overlapsAny(GRL1Beck_hg19_large, GR_MEInsCall))
mean(sapply(SampleColumns, function(x) sum(L1_1000G[idxFull,x] > 0))) / 
  mean(overlapsAny(GRL1Beck_hg19_large, GR_MEInsCall)) * 43 / 68

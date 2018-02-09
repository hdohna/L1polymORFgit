# The script below analyzes the locations of LINE-1 insertions in the human
# genome, distingushing between fragments, full-length and functional 
# insertions

#############################################
#                                           #
#    Source packages and set parameters     #
#                                           #
#############################################

# Load packages
library(ape)
library(seqinr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ShortRead)
library(csaw)

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.r')

# Path to L1 catalogue file (Created in script AddColumns2L1Catalog.R)
L1CataloguePath <- "D:/L1polymORF/Data/L1CatalogExtended.csv"

# Path to chain file for b37 to hg19 (obtained from ftp://gsapubftp-anonymous@ftp.broadinstitute.org/Liftover_Chain_Files/)
Chain_b37tohg19 <- "D:/L1polymORF/Data/b37tohg19.chain"

# Path to L1 GRanges from 1000 genome data
L1GRanges1000GenomesPath <- "D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData"

# Maximum L1 insertion length to be labeled a 'fragment'
MaxFragLength <- 5900

######################################
#                                    #
#    Read & process L1 catalog       #
#                                    #
######################################

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

# Create genomic ranges for catalog L1 in the reference
blnInRef <- (L1CatalogL1Mapped$end_HG38 - L1CatalogL1Mapped$start_HG38) > 6000 
L1CatalogGR_Ref     <- L1CatalogGR[blnInRef]
L1CatalogGR_nonRef  <- L1CatalogGR[!blnInRef]
blnZero <- L1CatalogL1Mapped$Activity == 0
L1CatalogGR_nonZero <- L1CatalogGR[!blnZero]

# Get fullength L1 that are not in the catalog
blnOverlapCatalog <- overlapsAny(L1RefGRFull, L1CatalogGR, minoverlap = 6000)
L1RefGRFullnotCat <- L1RefGRFull[!blnOverlapCatalog]
L1GRZero          <- c(L1RefGRFullnotCat, L1CatalogGR[blnZero])

# Create genomic ranges for catalog L1 in the reference
blnInRef_hg19   <- width(L1CatalogGR_hg19) > 6000 
L1CatalogGR_Ref_hg19 <- L1CatalogGR_hg19[blnInRef_hg19]

#######################################
#                                     #
#    Read & process recombination data  #
#                                     #
#######################################

cat("Processing 1000 Genome data\n")

# Read in file and create GRanges
# RecData <- read.delim("D:/L1polymORF/Data/hg19deCodeRecomb.txt")
RecData <- read.delim("D:/L1polymORF/Data/hg19RecombRate.txt")
Rec_GR <- makeGRangesFromDataFrame(RecData)
unique(seqnames(Rec_GR))
unique(seqnames(L1CatalogGR_Ref_hg19))

# Find overlaps to catalog L1 and create a vector of recombination values
RecL1CatOverlaps <- findOverlaps(L1CatalogGR_Ref_hg19, Rec_GR)
RecL1Cat         <- RecData$decodeAvg[RecL1CatOverlaps@to]  

# Analyze the effect of activity and recombination rate on frequency
any(duplicated(RecL1CatOverlaps@from))
L1CatalogRef_hg19 <- L1LiftoverList$L1CatalogWithHG19[blnInRef_hg19,]
LM <- lm(L1CatalogRef_hg19$Allele_frequency_Num[RecL1CatOverlaps@from] ~ 
           L1CatalogRef_hg19$ActivityNum[RecL1CatOverlaps@from] +
           RecL1Cat)
summary(LM)
cor.test(RecL1Cat, L1CatalogRef_hg19$ActivityNum[RecL1CatOverlaps@from])
plot(L1CatalogRef_hg19$ActivityNum[RecL1CatOverlaps@from], RecL1Cat)

# Find overlaps to fragment L1 and create a vector of recombination values
RecL1FragmOverlaps <- findOverlaps(L1FragmGR_hg19, Rec_GR)
RecL1Fragm <- RecData$decodeAvg[RecL1FragmOverlaps@to]  

# Test by sampling fragment L1s
NSamples <- 10000
SampledRecMeans <- sapply(1:NSamples, function(x){
  mean(sample(RecL1Fragm, length(RecL1Cat)))
})
sum(SampledRecMeans >= mean(RecL1Cat)) / NSamples



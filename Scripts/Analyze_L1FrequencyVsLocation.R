# The script below analyzes how the frequencies of L1s depends on their 
# locations with respect to genes and chromatin loops

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

#######################################
#                                     #
#    Read & process 1000 genome data  #
#                                     #
#######################################

cat("Processing 1000 Genome data\n")

# Load data
load(L1GRanges1000GenomesPath)

# Subset genomic ranges to get full-length and fragments in three different 
# frequency classes
blnFull1000G       <- L1_1000G_GR_hg19@elementMetadata@listData$InsLength >= 6000
GRL1Ins1000G_Full <- L1_1000G_GR_hg19[which(blnFull1000G)]

##################################
#                                #
#    Distance distributions      #
#                                #
##################################

##########
#  Defining auxiliary function
##########

# Auxiliary function to get distances to closest gene
Dist2ClosestGene <- function(GR){
  DistGeneObj <- distanceToNearest(GR, GRgenes, ignore.strand = T) 
  DistGeneObj@elementMetadata@listData$distance
}

# Auxiliary function to create genomic ranges for both sides of a loop
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


# Calculate distances from full-length and fragment L1 from 1000 genome data to
# closest gene
GRgenes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
L1DistGene_1000G <- Dist2Closest(L1_1000G_GR_hg19, GRgenes_hg19)

##########
#  Calculate distances to domain boundaries, catalog L1
##########

# List of files with domains and loops
LoopFiles <- list.files("D:/L1polymORF/Data/HiCData/", pattern = "looplist.txt",
                          full.names = T)
LoopFiles <- LoopFiles[! LoopFiles %in% grep(".gz", LoopFiles, value = T)]

# Get list of domain and loop ranges
LoopGRList     <- lapply(LoopFiles, getLoopGRs)
names(LoopGRList) <- sapply(LoopFiles, 
                              function(x) strsplit(x, "_")[[1]][2])
LoopDistList <- lapply(LoopGRList, function(x) {
  list(DistFragm = Dist2Closest(L1FragmGR_hg19, x$AllLoops),
       DistCatRef = Dist2Closest(L1CatalogGR_Ref_hg19, x$AllLoops),
       DistRefnotCat = Dist2Closest(L1RefGRFullnotCat_hg19, x$AllLoops))
})
LoopDist0 <- sapply(LoopDistList, function(x) {
  c(nrFragm = length(x$DistFragm),
    nrCatRef = length(x$DistCatRef),
    nrRefnotCat = length(x$DistRefnotCat),
    prop0Fragm = mean(x$DistFragm == 0),
    prop0CatRef = mean(x$DistCatRef == 0),
    prop0RefnotCat = mean(x$DistRefnotCat == 0),
    sum0Fragm = sum(x$DistFragm == 0),
    sum0CatRef = sum(x$DistCatRef == 0),
    sum0RefnotCat = sum(x$DistRefnotCat == 0))
})


# Summarize loop ranges (union or intersect)
LoopGR_Intersect <- LoopGRList[[1]]$AllLoops
for (i in 2:length(LoopGRList)){
  LoopGR_Intersect <- intersect(LoopGR_Intersect, LoopGRList[[i]]$AllLoops)
}

# Calculate distances from fragment L1, catalog L1 in the reference genome and
# full-length L1 that are not in the catalog
L1DistLoopIntersect_1000G  <- Dist2Closest(L1_1000G_GR_hg19, LoopGR_Intersect)

# Fit a model for frequency vs distance from genes and insertion length, 
# separately for fragments and full-length L1
blnFragm <- L1_1000G_reduced$InsLength < 6000
logFreq <- log(L1_1000G_reduced$Frequency)
logFreq[logFreq == -Inf] <- - 10

LM_fragm <- lm(logFreq[blnFragm] ~ 
                 L1_1000G_reduced$InsLength[blnFragm] +
                 L1DistGene_1000G[blnFragm])
summary(LM_fragm)

LM_full<- lm(logFreq[!blnFragm] ~ 
                 L1_1000G_reduced$InsLength[!blnFragm] +
                 L1DistGene_1000G[!blnFragm])
summary(LM_full)

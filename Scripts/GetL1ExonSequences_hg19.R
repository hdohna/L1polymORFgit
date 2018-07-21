##############################################
#
# General description:
#
#   The following script reads the repeat masker table, determines genomic
#   ranges, overlaps them with exons, and saves the L1 sequences overlapping
#   with exons as fasta file

# Input:
#
#     L1TableFileName: path to file that contains L1HS ranges in a table

# Output:
#   
#    : ...

##############################################

########################################
#                                      #
#  Source packages and set parameters  #
#                                      #
########################################

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.r')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Specify file paths
L1TableFileName   <- "D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv"
L1TableFileName   <- "D:/L1polymORF/Data/repeatsHg19_L1HS.csv"
L1TableFileName_hg38   <- "D:/L1polymORF_old/Data/repeatsHg38_L1HS.csv"
OutResults        <- 'D:/L1polymORF/Data/L1ExonSeqs.fastq'

#######################################
#                                     #
#    Read in data                     #
#                                     #
#######################################

# Read in table with L1 ranges
L1Table <- read.csv(L1TableFileName)

# Create GRanges objects with L1 Seqences
L1IRanges <- IRanges(start = L1Table$genoStart,
                     end = L1Table$genoEnd)
L1GRanges <- GRanges(seqnames = L1Table$genoName, ranges = L1IRanges,
                     strand = L1Table$strand)

# Get genomic ranges of genes
GeneGR <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
ExonGR <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
PromGR <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene)

GeneGR_L1OL  <- subsetByOverlaps(GeneGR, L1GRanges)
PromGR_L1OL  <- subsetByOverlaps(PromGR, L1GRanges)
writeLines(GeneGR_L1OL@elementMetadata@listData$gene_id, 
           con = "D:/L1polymORF/Data/L1IntersectGeneIDs_Refhg19")
writeLines(PromGR_L1OL@elementMetadata@listData$tx_id, 
           con = "D:/L1polymORF/Data/L1IntersectPromGeneIDs_Refhg19")


# Subset L1 ranges to get the ones o
L1exons <- subsetByOverlaps(L1GRanges, ExonGR, ignore.strand = T)
sum(overlapsAny(L1GRanges, exons(TxDb.Hsapiens.UCSC.hg19.knownGene)))

# Get L1 overlapping with genes
L1Table$blnL1GeneOL <- overlapsAny(L1GRanges, genes(TxDb.Hsapiens.UCSC.hg19.knownGene),
                                   ignore.strand = T)
L1Table$blnSameStrand <- NA
L1Table$blnSameStrand[L1GeneOL@from] <- 
  as.vector(strand(L1GRanges))[L1GeneOL@from] ==   
  as.vector(strand(GenesGR))[L1GeneOL@to]

# Smooth proportions of L1 overlapping with genes
L1Table$repStart
PropOLSmoothed_InsL <- supsmu(L1Table$repStart, 1*L1Table$blnL1GeneOL)
PropSameStrandSmoothed_InsL <- supsmu(L1Table$repStart, 1*L1Table$blnSameStrand)
plot(PropOLSmoothed_InsL$x, PropOLSmoothed_InsL$y, xlab = "L1 insertion start [bp]",
     ylab = "Proportion of L1 with positive selection signal", type = "l",
     col = "red", xlim = c(0, 500), ylim = c(0.2, 0.6))
lines(PropSameStrandSmoothed_InsL$x, PropSameStrandSmoothed_InsL$y, col = "blue")
points(L1Table$repStart[L1Table$blnL1GeneOL], rep(0.5, sum(L1Table$blnL1GeneOL)),
       col = "red")
points(L1Table$repStart[which(L1Table$blnL1GeneOL & L1Table$blnSameStrand)], 
       rep(0.5, sum(L1Table$blnL1GeneOL & L1Table$blnSameStrand, na.rm = T)),
       col = "blue")

################
# Analyze overlap between genes and L1 for hg19
################

# Get 5' UTRs 
UTR5P   <- fiveUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Get introns
introns <- intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene)
cat("Counting the number of genes per intron ... ")
NrIntrons <- sapply(introns, length)
cat("done!")
cat("Determinining ")
blnOLIntron <- overlapsAny(introns, L1GRanges)
idxOLIntron <- lapply(which(blnOLIntron), function(x) {
  idxOL <- which(overlapsAny(introns[[x]], L1GRanges))
  if (runValue(strand(introns[[x]])) == "-"){
    idxOL <- NrIntrons[x] - idxOL
  }
  idxOL
})



length(blnOLIntron)
overlapsAny(L1GRanges, unlist(UTR5P), ignore.strand = T)
L1Table$blnL1UTR5POL <- overlapsAny(L1GRanges, unlist(UTR5P), ignore.strand = T)
L1Table$blnL1intron  <- overlapsAny(L1GRanges, unlist(introns), ignore.strand = T)
mean(L1Table$blnL1UTR5POL)
mean(L1Table$blnL1GeneOL)

GenesGR  <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
L1GeneOL <- findOverlaps(L1GRanges, GenesGR, ignore.strand = T)

table(as.vector(strand(L1GRanges))[L1GeneOL@from], 
            as.vector(strand(GenesGR))[L1GeneOL@to])
fisher.test(as.vector(strand(L1GRanges))[L1GeneOL@from], 
            as.vector(strand(GenesGR))[L1GeneOL@to])

################
# Analyze overlap between genes and L1 for hg38
################

library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Read table with L1HS repeats
L1Table_hg38 <- read.csv(L1TableFileName_hg38)

# turn repeat table in genomic ranges
L1GRanges_hg38 <- makeGRangesFromDataFrame(L1Table_hg38, 
                                     seqnames.field = "genoName",
                                     start.field = "genoStart",
                                     end.field = "genoEnd")
GenesGR_hg38  <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
L1GeneOL_hg38 <- findOverlaps(L1GRanges_hg38, GenesGR_hg38)

table(as.vector(strand(L1GRanges_hg38))[L1GeneOL_hg38@from], 
      as.vector(strand(GenesGR_hg38))[L1GeneOL_hg38@to])
fisher.test(as.vector(strand(L1GRanges_hg38))[L1GeneOL_hg38@from], 
            as.vector(strand(GenesGR_hg38))[L1GeneOL_hg38@to])

L1Table_hg38$blnL1GeneOL <- overlapsAny(L1GRanges_hg38, GenesGR_hg38)
PropOLSmoothed_InsL <- supsmu(L1Table_hg38$repStart, 1*L1Table_hg38$blnL1GeneOL)
plot(PropOLSmoothed_InsL$x, PropOLSmoothed_InsL$y, xlab = "L1 insertion start [bp]",
     ylab = "Proportion of L1 within gene body", type = "l",
     col = "red", xlim = c(0, 6100))

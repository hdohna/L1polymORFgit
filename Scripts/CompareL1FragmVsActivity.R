# The script below compares the number of fragments mapped to a catalog element
# with L1 activity

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
library(ShortRead)
library(csaw)

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Path to L1 catalogue file 
L1CataloguePath <- "D:/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv"

# Path to bam file of mapped L1 fragments
BamFilePath <- 'D:/L1polymORF/Data/L1fragments_aln2catalog.sorted.bam'

# Path to L1 alignment file 
FastaFileName <- 'D:/L1polymORF/Data/L1CatalogueWithoutFlank_Sat_May_07_15-15-31_2016.fas'

# Minimum number of mapped fragments for L1 to be analyzed
MinFrags <- 20

# Maximum divergence for reads to be counted for L1 activity
MaxDiverge <- 0.002

######################################
#                                    #
#    Read & process data       #
#                                    #
######################################

# Read fasta file
L1Seq <- read.fasta(FastaFileName)
L1SeqLen <- sapply(L1Seq, length)

# Get reads from bam file aligned to L1HS concensus
L1CatBamGR <- GRanges(seqnames = names(L1Seq), 
                      IRanges(start = 1, end = L1SeqLen))
readParams  <- ScanBamParam(which = L1CatBamGR,
                            what = c('qname', "pos","cigar", "qwidth", "seq"), 
                            tag = "NM")
Reads       <- scanBam(BamFilePath, param = readParams)
names(Reads)
ReadsAccNrs <- sapply(names(Reads), function(x) strsplit(x, ":")[[1]][1])

# Get for each L1 element the proportion of nucs of mapped reads that differ 
# from reference L1
PropDiffList <- lapply(Reads, function(x) x$tag$NM / x$qwidth)
names(PropDiffList) <- ReadsAccNrs
hist(unlist(PropDiffList), breaks = seq(0, 0.2, 0.001))

# Get the number of reads
NrReads <- sapply(Reads, function(x) length(x$seq))
names(NrReads) <- ReadsAccNrs
sum(NrReads > 20)

# Get the number of close reads for each L1
NrCloseReads <- sapply(PropDiffList, function(x) sum(x <= MaxDiverge))

# Read in table with known L1 
L1Catalogue <- read.csv(L1CataloguePath, as.is = T)
AccMatch    <- match(L1Catalogue$Accession, names(NrReads))
L1Catalogue$NrReads <- NrReads[AccMatch]
L1Catalogue$NrCloseReads <- NrCloseReads[AccMatch]

# Get numeric activity
ActivityNum <- L1Catalogue$Activity
ActivityNum <- gsub("<", "", ActivityNum)
ActivityNum <- as.numeric(ActivityNum)
L1Catalogue$ActivityNum <- ActivityNum

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
L1CatalogL1Mapped <- L1Catalogue[blnL1Mapped & blnAllele1,]

# Create genomic ranges for catalog L1
L1CatalogGR <- GRanges(seqnames = L1CatalogL1Mapped$Chromosome,
   ranges = IRanges(start = pmin(L1CatalogL1Mapped$start_HG38,
                                 L1CatalogL1Mapped$end_HG38),
                    end = pmax(L1CatalogL1Mapped$start_HG38,
                               L1CatalogL1Mapped$end_HG38)),
                      strand = L1CatalogL1Mapped$Strand)

# Plot activity vs number reads
plot(L1Catalogue$Activity, L1Catalogue$NrReads)
plot(L1Catalogue$ActivityNum * as.numeric(L1Catalogue$Allele_frequency),
     L1Catalogue$NrReads)
plot(L1Catalogue$ActivityNum * as.numeric(L1Catalogue$Allele_frequency),
     L1Catalogue$NrCloseReads)

# Create distance classes
DiffClasses <- seq(0, 0.2, 0.002)
idxSomeFrags <- which(NrReads >= MinFrags)
Cols <- rainbow(length(idxSomeFrags))
plot(DiffClasses, DiffClasses, ylim = c(0, 40), type = "n",
     xlim = c(0, 0.05), xlab = "Difference to full-length",
     ylab = "Count")

# Create a matrix that counts for each full-length L1 that has some fragments 
# associated the number of associated fragments per distance class
CountMat <- sapply(1:length(idxSomeFrags), function(i){
  x <- PropDiffList[[idxSomeFrags[i]]]
  H <- hist(x, breaks = DiffClasses, plot = F)
  lines(H$mids, H$counts, col = Cols[i], type = "s")
  H$counts
})
colnames(CountMat) <- names(NrReads)[idxSomeFrags]
rownames(CountMat) <- DiffClasses[-1]


# Get ranges of genes and L1 overlapping with genes
GRgenes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
L1OverlapGene <- overlapsAny(L1CatalogGR, GRgenes)

# Calculate product of activity and frequency and plot it against nr mapped
# reads
ActFreq <- L1CatalogL1Mapped$ActivityNum * 
   as.numeric(L1CatalogL1Mapped$Allele_frequency)
plot(ActFreq, L1CatalogL1Mapped$NrCloseReads)
points(ActFreq[L1OverlapGene], L1CatalogL1Mapped$NrCloseReads[L1OverlapGene],
       col = "red")
plot(L1CatalogL1Mapped$ActivityNum, L1CatalogL1Mapped$NrReads)
points(L1CatalogL1Mapped$ActivityNum[L1OverlapGene], 
       L1CatalogL1Mapped$NrReads[L1OverlapGene],
       col = "red")


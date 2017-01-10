# The script below reads a list of reads aligned to L1, a coverage and quantile 
# matrix (created in script 'CalcCoverMatReadList.R') and creates new L1 
# reference with known and new L1 stumps

##############
# Source prerequisites
##############

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(ShortRead)
library(Rsamtools)
library(rtracklayer)
library(csaw)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

##############
# Set parameters
##############

# Set parameters for range to determine start and end of 5' region of L1 stumps
FivePStart <- 50
FivePEnd   <- 2000

# Set file paths
CoveragePlotPath      <- 'D:/L1polymORF/Figures/L1InsertionCoverage_NA12878_PacBio.pdf'
CoverDataPath         <- 'D:/L1polymORF/Data/L1_NA12878_PacBio_Coverage.RData'
L1_1000GenomeDataPath <- "D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData"
OutFolderName_NonRef  <- "D:/L1polymORF/Data/BZ_NonRef"
NewL1RefOutPath       <- "D:/L1polymORF/Data/L1RefPacBioNA12878.fa"

##############
# Load data
##############

# Load ranges
load("D:/L1polymORF/Data/L1RefRanges_hg19.RData")
load("D:/L1polymORF/Data/BZ_L1Ranges.RData")
load(L1_1000GenomeDataPath)
load(CoverDataPath)

# Read in table with known L1 
L1Catalogue <- read.csv("D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv", as.is = T)

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
L1CatalogL1Mapped <- L1Catalogue[which(blnL1Mapped & blnAllele1),]

# Lift coordinates and get genomic ranges for catalog L1 on hg19
LiftOverList <- LiftoverL1Catalog(L1CatalogL1Mapped,
                                  ChainFilePath = "D:/L1polymORF/Data/hg38ToHg19.over.chain")
L1CatalogGR <- LiftOverList$GRCatalogue_hg19# Specify folder 

# Read in L1 consensus sequence
L1ConsSeq <- read.fasta("D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fas")
L1ConsSeq <- toupper(L1ConsSeq[[1]])

# Read repeat table and subset to get only L1HS rows with fragment size below 
# MaxFragLength
RepeatTable <- read.csv("D:/L1polymORF/Data/repeatsHg38_L1HS.csv")
RepeatTable <- read.delim("D:/L1polymORF/Data/repeatL1hg19")
RepeatTable <- RepeatTable[RepeatTable$repName == "L1HS", ]

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
L1RefGR <- GRanges(seqnames = RepeatTable$genoName,
                   ranges = IRanges(start = RepeatTable$genoStart,
                                    end = RepeatTable$genoEnd),
                   strand = RepeatTable$strand)
L1RefGRFull <- L1RefGR[width(L1RefGR) > 6000]

L1RefSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, L1RefGRFull)
L1RefSeq <- as.character(L1RefSeq)
StrandV  <- as.vector(strand(L1RefGRFull))
L1RefSeqMat <- sapply(1:length(L1RefSeq), function(i) {
  L1V <- strsplit(L1RefSeq[i], "")[[1]]
  if (StrandV[i] == "+"){
    L1V[1:FivePEnd]
  } else {
    L1V[(length(L1V) - FivePEnd + 1):length(L1V)]
  }
})
dim(L1RefSeqMat)
colnames(L1RefSeqMat) <- paste(as.vector(seqnames(L1RefGRFull)), start(L1RefGRFull),
                               end(L1RefGRFull), "Ref", sep = "_")

# get names of newly created bam files
FileNames <- list.files(OutFolderName_NonRef, pattern = ".bam",
                        full.names = T)
FileNames <- FileNames[-grep(".bam.", FileNames)]

# Collect information on insertion that fullfill a certain minimum criterion
idx5P <- which(sapply(1:nrow(CoverMat), function(x) all(CoverMat[x, FivePStart:FivePEnd] > 0)))
idxFull <- which(sapply(1:nrow(CoverMat), function(x) all(CoverMat[x, FivePStart:6040] > 0)))
FullL1Info <- t(sapply(FilesWithReads[idx5P], function(x){
  FPathSplit <- strsplit(x, "/")[[1]]
  FName      <- FPathSplit[length(FPathSplit)]
  FName      <- substr(FName, 1, nchar(FName) - 4)
  strsplit(FName, "_")[[1]]
}))
FullL1Info <- base::as.data.frame(FullL1Info, stringsAsFactors = F)
colnames(FullL1Info) <- c("chromosome", "idx")
FullL1Info$Max3P <- sapply(idx5P, function(x) max(which(CoverMat[x, ] > 0)))
FullL1Info$idx   <- as.numeric(FullL1Info$idx)
FullL1Info$start <- start(IslGRanges_reduced[FullL1Info$idx])
FullL1Info$end   <- end(IslGRanges_reduced[FullL1Info$idx])
FullL1Info$cover   <- CoverMat[idx5P, 100]
FilesWithReads[idxFull]

# Create genomic ranges for L1 insertions and compare with known L1 ranges
GRL1Capture <- makeGRangesFromDataFrame(FullL1Info)
GRL1Capture100 <- resize(GRL1Capture, 100, fix = "center")
L1CatalogGR100 <- resize(L1CatalogGR, 100, fix = "center")
GRL1Ins1000G100 <- resize(GRL1Ins1000G, 100, fix = "center")
if (any(c(sum(overlapsAny(GRL1Capture100, L1CatalogGR100)),
          sum(overlapsAny(GRL1Capture100, GRL1Ins1000G100)),
          sum(overlapsAny(IslGRanges_reduced, L1CatalogGR100)),
          sum(overlapsAny(IslGRanges_reduced, GRL1Ins1000G100))) > 0)){
  cat("Some newly found L1 insertion overlap with catalog or 1000 Genome elements!\n")
}

# Get sequences of a particular insertion
NewL1Seq <- sapply(1:length(idx5P), function(i){
  x <- idx5P[i]
  cat("Processing insertion", i, "of", length(idx5P), "\n")
  RL <- ReadListPerPeak[[x]]
  idxReads <- which(RL$pos < FivePEnd)
  
  # Loop through reads and create a matrix of aligned reads
  ReadMat  <- sapply(idxReads, function(i) {
    SeqV   <- SeqFromCigar(RL$cigar[i], RL$seq[i])
    SeqEnd  <- min(FivePEnd - RL$pos[i] + 1, length(SeqV))
    Prepend <- rep("-", RL$pos[i] - 1)
    NrAppend <- FivePEnd - RL$pos[i] + 1 - SeqEnd
    Append  <- rep("-", NrAppend)
    SV <- c(Prepend, SeqV[1:SeqEnd], Append)
  })
  
  # Create a consensus sequence
  ConsensSeq  <- c()
  ConsensProp <- c()
  for(x in 1:nrow(ReadMat)){
    NucCount    <- table(ReadMat[x,])
    ConsensSeq  <- c(ConsensSeq,names(NucCount)[which.max(NucCount)])
    ConsensProp <- c(ConsensProp,max(NucCount)/sum(NucCount))
  }
  
  # Replace nucleotides that are variable 
  blnReplace <- (ConsensSeq == "-") | ConsensProp < 0.6
  ConsensSeq[blnReplace] <- L1ConsSeq[blnReplace]
  ConsensSeq
})
colnames(NewL1Seq) <- paste(FullL1Info$chromosome, FullL1Info$start,
                            FullL1Info$end, "New", sep = "_")

# Calculate the number of nucleaotide diefferences between different L1 stumps
DiffMat <- sapply(1:ncol(NewL1Seq), function(x) {
  sapply(1:ncol(NewL1Seq), function(y){
     sum(NewL1Seq[ , x] != NewL1Seq[ , y])
  })
})

# Count differences to consensus sequence
Diff2L1Consens <- sapply(1:ncol(NewL1Seq), function(x) {
    sum(NewL1Seq[ , x] != L1ConsSeq[1:FivePEnd])
})
hist(Diff2L1Consens, breaks = seq(-5, 105, 5))
sum(Diff2L1Consens == 0)

# Determine indices of sequences that are identical with another
diag(DiffMat) <- NA
max(DiffMat, na.rm = T)
min(DiffMat, na.rm = T)
sum(DiffMat == 0, na.rm = T)
idxIdentical <- which(DiffMat == 0, arr.ind = T)
idxIdentical <- unique(as.vector(idxIdentical))

# Get L1 5' sequences that are non-reference
L1Catalogue <- LiftOverList$L1CatalogWithHG19
blnNonRef   <- (L1Catalogue$end_HG19 - L1Catalogue$start_HG19) <= 5000
idxNonRef   <- which(blnNonRef)
L1SeqNonRef <- L1Catalogue$L1Seq[idxNonRef]

# Turn all L1 seq 
# Pattern <- paste(L1ConsSeq[100:200], collapse = "")
# pMatch <- vmatchPattern(Pattern, L1SeqNonRef, max.mismatch = 5)
# sapply(pMatch, length)
L1NonRef <- sapply(L1SeqNonRef, function(x) strsplit(x, "")[[1]][1:FivePEnd])
colnames(L1NonRef) <- paste(L1Catalogue$Chromosome[idxNonRef], 
                            L1Catalogue$start_HG19[idxNonRef],
                            L1Catalogue$end_HG19[idxNonRef], 
                            L1Catalogue$Accession[idxNonRef], sep = "_")


# Put all the different L1 sequences together
L1StumpRef <- cbind(NewL1Seq[,-idxIdentical[-1]], L1NonRef, L1RefSeqMat)
L1StumpRef <- t(L1StumpRef)
dim(L1StumpRef)
#write.dna(L1StumpRef, NewL1RefOutPath, format = "fasta")
L1StumpList <- lapply(1:nrow(L1StumpRef), function(x) L1StumpRef[x,])
write.fasta(L1StumpList, rownames(L1StumpRef), NewL1RefOutPath)


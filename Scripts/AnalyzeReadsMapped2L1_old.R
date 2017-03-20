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

############################
#                          #
#      Set parameters      #
#                          #
############################

# Set parameters for range to determine start and end of 5' region of L1 stumps
# and the minimum coverage in that region
FivePStart <- 50
FivePEnd   <- 2000
MinCover   <- 1

# Set parameter to determine whether a L1 insertion is potentially full-length
# (it is not full-length if it has reads that are clipped by MinClip or more,
# at position MaxL1Pos or less)
MinClip  <- 500
MaxL1Pos <- 5500

# Set parameter for minimum proportion of a polymorphism among reads to be called
MinPolyProp <- 0.6

# Specify flank size for sequences flanking L1 to be included in alignment output
FlankSize <- 500

# Minimum number of bp clipped of a read for it to be included in alignment
MinClipAlign <- 1000

# Set file paths
CoveragePlotPath      <- 'D:/L1polymORF/Figures/L1InsertionCoverage_NA12878_PacBio.pdf'
CoverDataPath         <- 'D:/L1polymORF/Data/L1_NA12878_PacBio_Coverage.RData'
L1_1000GenomeDataPath <- "D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData"
OutFolderName_NonRef  <- "D:/L1polymORF/Data/BZ_NonRef"
NewL1RefOutPath       <- "D:/L1polymORF/Data/L1RefPacBioNA12878_DelRemoved.fa"
GenomeBamPath         <- "D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_subreads_hg19masked.sorted.bam"
ResultPath            <- "D:/L1polymORF/Data/ReadsMapped2L1Info.RData"
FastaFilePath         <- "D:/L1polymORF/Data/"

################################
#                              #
#  Load and preprocess data    #
#                              #
################################

######
# Load ranges
######

# Load ranges
load("D:/L1polymORF/Data/L1RefRanges_hg19.RData")
load("D:/L1polymORF/Data/BZ_L1Ranges.RData")
load(L1_1000GenomeDataPath)
load(CoverDataPath)

######
# Read and process L1 catalog
######

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

######
# Read and process L1 consensus sequence
######

# Read in L1 consensus sequence
L1ConsSeq <- read.fasta("D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fas")
L1ConsSeq <- toupper(L1ConsSeq[[1]])
L1ConsSeqDNAst <- DNAString(paste(L1ConsSeq, collapse = ""))

######
# Read and process repeatmasker table
######

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

###############################################
#                                             #
#  Collect info on promising L1 insertions    #
#            (FullL1Info)                     #
#                                             #
###############################################

# get names of newly created bam files
FileNames <- list.files(OutFolderName_NonRef, pattern = ".bam",
                        full.names = T)
FileNames <- FileNames[-grep(".bam.", FileNames)]

# Collect information on insertion that fullfill a certain minimum criterion
idx5P <- which(sapply(1:nrow(CoverMat), function(x) {
  all(CoverMat[x, FivePStart:FivePEnd] >= MinCover)
}))
cat("******  ", length(idx5P), "insertions have a minimum coverage of", 
MinCover, "in 5' region from", FivePStart, "to", FivePEnd, "   *********\n")
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
idxFull2 <- which(rownames(FullL1Info) %in% FilesWithReads[idxFull])

# Get indices of insertions that are not full-length
FullL1Info$PotentialFullLength <- sapply(1:length(idx5P), function(i){
  x <- idx5P[i]
  RL <- ReadListPerPeak[[x]]
  primMap <- RL$flag <= 2047
  RL <- lapply(RL, function(y) y[primMap])
  ReadLenghsL1 <- sapply(RL$cigar, ReadLengthFromCigar)
  LRClipped <- sapply(RL$cigar, NrClippedFromCigar)
  !any((ReadLenghsL1 <= MaxL1Pos) & (LRClipped[2,] > MinClip))
})
FullL1Info$PotentialFullLength[idxFull2] <- T
cat("******  ", sum(FullL1Info$PotentialFullLength, na.rm = T), 
    "insertions are potentially full length *********\n")

# Set up insertion descriptors
FullL1GR <- makeGRangesFromDataFrame(FullL1Info)
FullL1Info$L1InsertionPosition.median <- NA
FullL1Info$L1InsertionPosition.min    <- NA
FullL1Info$L1InsertionPosition.max    <- NA
FullL1Info$L1Strand                   <- NA
FullL1Info$L15PTransdSeq.median       <- NA
FullL1Info$L15PTransdSeq.min          <- NA
FullL1Info$L15PTransdSeq.max          <- NA
FullL1Info$NrSupportReads             <- NA

# Loop through potential full-length insertions
for (i in which(FullL1Info$PotentialFullLength)){
  x <- idx5P[i]
  
  # Get reads mapped to L1, retain only primary reads and get the number of bp
  # clipped on the left
  RL <- ReadListPerPeak[[x]]
  primMap     <- RL$flag <= 2047
  RL          <- lapply(RL, function(y) y[primMap])
  LRClippedL1 <- sapply(RL$cigar, NrClippedFromCigar)
  L1Length    <- width(RL$seq) - LRClippedL1[1,]
  
  # Get reads mapped to the genome on current locus
  ScanParam <- ScanBamParam(what = scanBamWhat(), which = FullL1GR[i])
  GenomeRL  <- scanBam(GenomeBamPath,  param = ScanParam)
  LRClippedG <- sapply(RL$cigar, NrClippedFromCigar)
  ReadMatch <- match(RL$qname, GenomeRL[[1]]$qname)
  
  # Loop through reads an collect information on number of bps clipped on
  # reads mapped to genome, the insertion location and the insertion strand
  L1Info <- sapply(1:length(ReadMatch), function(y) {
    idxGR <- ReadMatch[y]
    LRClipped <- NrClippedFromCigar(GenomeRL[[1]]$cigar[idxGR])
    if (LRClipped[1] > LRClipped[2]) {
      bpClipped <- LRClipped[1]
      L1Strand  <- 0 # 0 corresponds to negative strand
      InsPos    <- GenomeRL[[1]]$pos[idxGR]
    } else {
      bpClipped <- LRClipped[2]
      L1Strand  <- 1 # 1 corresponds to positive strand
      InsPos <- GenomeRL[[1]]$pos[idxGR] + 
        ReadLengthFromCigar(GenomeRL[[1]]$cigar[idxGR])
    }
    c(bpClipped = bpClipped, L1Strand = L1Strand, InsPos = InsPos)
  })
  
  # Fill in info about L1 insertion in FullL1Info
  if (mean(L1Info["L1Strand",]) > 0.5) {
    FullL1Info$L1Strand[i] <- "+"
  } else {
    FullL1Info$L1Strand[i] <- "-"
  }
  FullL1Info$L1InsertionPosition.median[i] <- median(L1Info["InsPos",])
  FullL1Info$L1InsertionPosition.min[i]    <- min(L1Info["InsPos",])
  FullL1Info$L1InsertionPosition.max[i]    <- max(L1Info["InsPos",])
  
  # Get length of transduced sequence and number of supporting reads
  LengthTransduced <- -(L1Length - L1Info["bpClipped",])
  FullL1Info$L15PTransdSeq.median[i] <- median(LengthTransduced)
  FullL1Info$L15PTransdSeq.min[i]    <- min(LengthTransduced)
  FullL1Info$L15PTransdSeq.max[i]    <- max(LengthTransduced)
  FullL1Info$NrSupportReads[i]       <- length(ReadMatch)
}

# Loop through potential full-length insertions and 
rownames(FullL1Info)[which(FullL1Info$Max3P >= 6020)]
idxFull2 <- which(rownames(FullL1Info) %in% FilesWithReads[idxFull])
i <- 68
for (i in idxFull2){

  # Create a file name for output fasta file
  FastaFile <- paste(FastaFilePath, FullL1Info$chromosome[i], "_",FullL1Info$idx[i],
                     ".fas", sep = "")
  
  # Get insertion location and create a reference sequence with
  InsLoc     <- FullL1Info$L1InsertionPosition.median[i]
  L1Str      <- FullL1Info$L1Strand[i]
  L1ConsSeqLocal <- L1ConsSeqDNAst
  if (L1Str == "-") L1ConsSeqLocal <- reverseComplement(L1ConsSeqLocal)
  TransDL    <- FullL1Info$L15PTransdSeq.median[i]
  LeftFlangR <- GRanges(FullL1Info$chromosome[i], 
                        IRanges(start = InsLoc - FlankSize + 1,
                                   end = InsLoc + FlankSize))
  RightFlangR <- GRanges(FullL1Info$chromosome[i], IRanges(start = InsLoc + 1,
                            end = InsLoc + FlankSize))
  LeftFlankSeq  <- getSeq(BSgenome.Hsapiens.UCSC.hg19, LeftFlangR)
  RightFlankSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, RightFlangR)
  RefSeq <- c(LeftFlankSeq[[1]], L1ConsSeqLocal, RightFlankSeq[[1]])
  RefSeqCharV <- s2c(as.character(RefSeq))
  
  # Get reads mapped to the genome on current locus
  ScanParam  <- ScanBamParam(what = scanBamWhat(), which = FullL1GR[i])
  GenomeRL   <- scanBam(GenomeBamPath,  param = ScanParam)
  LRClippedG <- sapply(GenomeRL[[1]]$cigar, NrClippedFromCigar)
  idxEnoughClipped <- which(
    pmax(LRClippedG[1,], LRClippedG[2,]) > MinClipAlign)

  # Loop through reads an collect sequence
  SeqList <- vector(mode = "list", length = length(idxEnoughClipped) + 1)
  SeqList[[1]] <- RefSeqCharV
  idxL <- 1
  for (y in idxEnoughClipped) {
    LRClipped  <- LRClippedG[,y]
    Seq        <- GenomeRL[[1]]$seq[y]
    if (LRClipped[1] > LRClipped[2]) {
      SeqEnd   <- min(width(Seq), LRClipped[1] + TransDL + FlankSize)
      Seq2Add  <- subseq(Seq, 1, SeqEnd)
    } else {
      SeqStart <- max(width(Seq), width(Seq) - LRClipped[2] - TransDL - FlankSize)
      Seq2Add  <- subseq(Seq, SeqStart, width(Seq))
    }
#    if (GenomeRL[[1]]$strand[y] == "-") Seq2Add <- reverseComplement(Seq2Add)
    idxL <- idxL + 1
    SeqList[[idxL]] <- s2c(as.character(Seq2Add))
  }
  names(SeqList) <- c("RefSeq", GenomeRL[[1]]$qname[idxEnoughClipped])
  names(SeqList) <- gsub("/", "_", names(SeqList))
  write.fasta(SeqList, names = names(SeqList), file.out = FastaFile)
  
  AlignedFile <- gsub(".fas", "_aligned.fas", FastaFile)
  cat("Aligning", FastaFile, "\n")
  run_MUSCLE(InputPath = FastaFile, OutputPath = AlignedFile) 
}

# Create genomic ranges for L1 insertions and compare with known L1 ranges
GRL1Capture <- makeGRangesFromDataFrame(FullL1Info)
GRL1Capture100 <- resize(GRL1Capture, 100, fix = "center")
L1CatalogGR100 <- resize(L1CatalogGR, 100, fix = "center")
GRL1Ins1000G100 <- resize(GRL1Ins1000G, 100, fix = "center")
if (any(c(sum(overlapsAny(GRL1Capture100, L1CatalogGR100)),
          sum(overlapsAny(GRL1Capture100, GRL1Ins1000G100))) > 0)){
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
  ConsensSeq  <- L1ConsSeq[1:FivePEnd]
  ConsensProp <- rep(1, FivePEnd)
  AllDel <- sapply(1:nrow(ReadMat), function(x) all(ReadMat[x,] == "-"))
  ConsensSeq[AllDel] <- "-"
  for(x in which(!AllDel)){
    NucCount    <- table(ReadMat[x,])
    NucCount    <- NucCount[names(NucCount) != "-"]
    ConsensSeq[x]  <- names(NucCount)[which.max(NucCount)]
    ConsensProp[x] <- max(NucCount)/sum(NucCount)
  }
  
  # Replace nucleotides that are variable 
  idxReplace <- which(ConsensProp < MinPolyProp)
  ConsensSeq[idxReplace] <- L1ConsSeq[idxReplace]
  ConsensSeq
})
colnames(NewL1Seq) <- paste(FullL1Info$chromosome, FullL1Info$start,
                            FullL1Info$end, "New", sep = "_")

# Calculate the number of nucleotide diefferences between different L1 stumps
DiffMat <- sapply(1:ncol(NewL1Seq), function(x) {
  cat("Calculating differences for sequence", x, "of", ncol(NewL1Seq), "\n")
  sapply(1:ncol(NewL1Seq), function(y){
     sum(NewL1Seq[ , x] != NewL1Seq[ , y])
  })
})

# Determine indices of sequences that are identical with another
diag(DiffMat) <- NA
max(DiffMat, na.rm = T)
min(DiffMat, na.rm = T)
mean(DiffMat, na.rm = T)
sum(DiffMat == 0, na.rm = T)
idxIdentical <- which(DiffMat == 0, arr.ind = T)
idxIdentical <- unique(as.vector(idxIdentical))

# Count differences to consensus sequence
Diff2L1Consens <- sapply(1:ncol(NewL1Seq), function(x) {
    sum(NewL1Seq[ , x] != L1ConsSeq[1:FivePEnd])
})
hist(Diff2L1Consens, breaks = seq(-5, 1005, 5))
sum(Diff2L1Consens == 0)
mean(Diff2L1Consens)
NrDel <- sapply(1:ncol(NewL1Seq), function(x) {
  sum(NewL1Seq[ , x] == "-")
})
rownames(FullL1Info)[which.max(NrDel)]

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


# Put all the different L1 sequences together, remove consistent deletions and
# save as fasta file
L1StumpRef <- cbind(NewL1Seq, L1NonRef, L1RefSeqMat)
L1StumpList <- lapply(1:ncol(L1StumpRef), function(x) {
  c(L1StumpRef[L1StumpRef[,x] != "-",x], L1ConsSeq[(FivePEnd + 1):length(L1ConsSeq)])
})
cat("*****   Saving new L1 references as file", NewL1RefOutPath, "  *****\n")
write.fasta(L1StumpList, colnames(L1StumpRef), NewL1RefOutPath)

# Save info about insertions
cat("*****   Saving image to", ResultPath, "  *****\n")
save.image(ResultPath)


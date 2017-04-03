# The script below reads a list of reads aligned to L1, a coverage and quantile 
# matrix (created in script 'CalcCoverMatReadList.R') and creates new L1 
# reference with known and new L1 stumps

##############
# Source prerequisites
##############

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(smooth)
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

# Width of a motif to determine which side of a clipped read maps to L1
MotifWidth  <- 50
MotifOffSet <- 60

# Specify flank size for sequences flanking L1 to be included in alignment output
FlankSize <- 2000

# Minimum number of bp clipped of a read for it to be included in alignment
MinClipAlign <- 1000

# Minimum length of full-length L1
MinFullL1 <- 6020

# Set file paths
CoveragePlotPath      <- 'D:/L1polymORF/Figures/L1InsertionCoverage_NA12878_PacBio.pdf'
CoverDataPath         <- 'D:/L1polymORF/Data/L1_NA12878_PacBio_Coverage.RData'
L1_1000GenomeDataPath <- "D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData"
OutFolderName_NonRef  <- "D:/L1polymORF/Data/BZ_NonRef"
NewL1RefOutPath       <- "D:/L1polymORF/Data/L1RefPacBioNA12878_DelRemoved.fa"
GenomeBamPath         <- "D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_subreads_hg19masked.sorted.bam"
#GenomeBamPath         <- "D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_subreads_hg19withL1.sorted.bam"
ResultPath            <- "D:/L1polymORF/Data/ReadsMapped2L1Info.RData"
FastaFilePath         <- "D:/L1polymORF/Data/"
L1InsertionFastaPath  <- paste(FastaFilePath, "L1InsertionWithFlank", FlankSize,
                               "bp.fas", sep = "")

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

# Turn into a DNAstring and create reverse complement
L1ConsSeqDNAst    <- DNAString(paste(L1ConsSeq, collapse = ""))
L1ConsSeqDNAst_RV <- reverseComplement(L1ConsSeqDNAst)
L1CharV           <- s2c(as.character(L1ConsSeqDNAst))
L1CharV_RV        <- s2c(as.character(L1ConsSeqDNAst_RV))

######
# Read and process repeatmasker table
######

# Read repeat table and subset to get only L1HS rows with fragment size below 
# MaxFragLength
RepeatTable <- read.csv("D:/L1polymORF/Data/repeatsHg38_L1HS.csv")
#RepeatTable <- read.delim("D:/L1polymORF/Data/repeatL1hg19")
RepeatTable <- read.delim("D:/L1polymORF/Data/repeatsHg19_L1HS")
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

idxFull2 <- which(rownames(FullL1Info) %in% FilesWithReads[idxFull])

# Determine whether L1 might be potentially full-length
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
FullL1Info$L1StrandNum                <- NA
FullL1Info$L1Strand                   <- NA
FullL1Info$L1Start.median             <- NA
FullL1Info$L1Start.min                <- NA
FullL1Info$L1Start.max                <- NA
FullL1Info$L1End.median               <- NA
FullL1Info$L1End.min                  <- NA
FullL1Info$L1End.max                  <- NA
FullL1Info$L15PTransdSeq.median       <- NA
FullL1Info$L15PTransdSeq.min          <- NA
FullL1Info$L15PTransdSeq.max          <- NA
FullL1Info$L13PTransdSeq.median       <- NA
FullL1Info$L13PTransdSeq.min          <- NA
FullL1Info$L13PTransdSeq.max          <- NA
FullL1Info$NrSupportReads             <- NA
FullL1Info$NrReadsCover5P             <- NA
FullL1Info$NrReadsCover3P             <- NA
FullL1Info$NrReadsReach5P             <- NA
FullL1Info$NrReadsReach3P             <- NA

# Loop through potential full-length insertions and create a list of 
# read-specific information
idxPotentialFull <- which(FullL1Info$PotentialFullLength)
idxPotentialFull <- 1:nrow(FullL1Info)

# Generate a list of matched reads mapped to genome and to L1
MatchedReadList <- lapply(idxPotentialFull, function(i) {
  
  # Get position of current entry in ReadListPerPeak
  x <- idx5P[i]
  
  # Get reads mapped to L1, retain only primary reads and get the number of bp
  # clipped on the left
  RL          <- ReadListPerPeak[[x]]
  primMap     <- RL$flag <= 2047 & !(is.na(RL$pos))
  RL          <- lapply(RL, function(y) y[primMap])
  LRClippedL1 <- sapply(RL$cigar, NrClippedFromCigar, USE.NAMES = F)
  L1Length    <- width(RL$seq) - LRClippedL1[1,]
  
  # Get reads mapped to the genome on current locus
  ScanParam <- ScanBamParam(what = scanBamWhat(), which = FullL1GR[i])
  GenomeRL  <- scanBam(GenomeBamPath,  param = ScanParam)
  primMap     <- GenomeRL[[1]]$flag <= 2047 & !(is.na(GenomeRL[[1]]$pos))
  GenomeRL[[1]]    <- lapply(GenomeRL[[1]], function(y) y[primMap])
  LRClippedG <- sapply(RL$cigar, NrClippedFromCigar)
  ReadMatch <- match(RL$qname, GenomeRL[[1]]$qname)
  if (any(is.na(ReadMatch))){
    browser()
  }
  
  # Generate a list of reads mapped to genome, matching the list of reads 
  # mapped to L1
  GenomeMatchRL <- lapply(GenomeRL[[1]], function(y) y[ReadMatch])
  list(L1ReadList = RL, GenomeReadList = GenomeMatchRL)
})

# Loop through list of matching reads and collect info
ReadInfoList <- lapply(MatchedReadList, function(RLL) {

  # Determine whether reads are on the same strand  
  blnSameStrand <- RLL$GenomeReadList$strand == RLL$L1ReadList$strand
  
  # Loop through reads an collect information on number of bps clipped on
  # reads mapped to genome, the insertion location and the insertion strand
  L1Info <- t(sapply(1:length(RLL$GenomeReadList$pos), function(y) {
    
    # Get sequence and nr clipped of read mapped to L1 
    L1Seq        <- RLL$L1ReadList$seq[y]
    LRClipped_L1 <- NrClippedFromCigar(RLL$L1ReadList$cigar[y])
    
    # Get start and end position on L1
    L1Start <- RLL$L1ReadList$pos[y]
    L1End   <- L1Start + ReadLengthFromCigar(RLL$L1ReadList$cigar[y])
    
    # Get sequence and nr clipped of read mapped to genome 
    GSeq        <- RLL$GenomeReadList$seq[y]
    LRClipped_G <- NrClippedFromCigar(RLL$GenomeReadList$cigar[y])
    
    # Determine the clipped sequence of read mapped to genome 
    LeftClippedSeq  <- subseq(GSeq, 1, LRClipped_G[1])
    RightClippedSeq <- subseq(GSeq, width(GSeq) - LRClipped_G[2], width(GSeq))
    
    # Get a part of read mapped to L1 and determine whether it is in the left
    # or right clipped sequence of read mapped to genome
    L1Motif         <- subseq(L1Seq, LRClipped_L1[1] + MotifOffSet,
                              LRClipped_L1[1] + MotifOffSet + MotifWidth)
    if (!blnSameStrand[y]) {L1Motif <- reverseComplement(L1Motif)}
    motifMatchLeft  <- matchPattern(L1Motif[[1]], LeftClippedSeq[[1]])
    motifMatchRight <- matchPattern(L1Motif[[1]], RightClippedSeq[[1]])
    blnL1OnLeft  <- length(motifMatchLeft)  > 0
    blnL1OnRight <- length(motifMatchRight) > 0
    
    # Identify position of L1 insertion
    InsPos <- NA
    if (blnL1OnLeft){
      InsPos <- RLL$GenomeReadList$pos[y]
    } 
    if (blnL1OnRight){
      InsPos <- RLL$GenomeReadList$pos[y] + 
        ReadLengthFromCigar(RLL$GenomeReadList$cigar[y])
    } 
    if (blnL1OnLeft & blnL1OnRight){
      warning("L1 motif matches left and right clipped sequence in", i, "\n", immediate. = T)
    }
    if ((length(motifMatchLeft) == 0) & (length(motifMatchRight) == 0)){
       browser()
       warning("L1 motif matches no clipped sequence in", i, "\n")
    }
    # Scenario: 0 = L1 on negative strand and read maps on left side of L1
    #           1 = L1 on positive strand and read maps on left side of L1
    #           2 = L1 on negative strand and read maps on right side of L1
    #           3 = L1 on positive strand and read maps on right side of L1
    Scenario <- blnSameStrand[y] + 2 * blnL1OnLeft
    
    # Get start and stop indices of transduced sequence of 5' and 3' end
    TD5Pstart <- switch(as.character(Scenario), 
                        '0' = width(GSeq) - LRClipped_L1[1], 
                        '1' = width(GSeq) - LRClipped_G[2],
                        '2' = width(GSeq) - LRClipped_L1[1],
                        '3' = 1)
    TD5Pend   <- switch(as.character(Scenario), 
                        '0' = width(GSeq), 
                        '1' = LRClipped_L1[1],
                        '2' = LRClipped_G[1],
                        '3' = LRClipped_L1[1])
    TD3Pstart <- switch(as.character(Scenario), 
                        '0' = width(GSeq) - LRClipped_G[2], 
                        '1' = width(GSeq) - LRClipped_L1[2],
                        '2' = 1,
                        '3' = width(GSeq) - LRClipped_L1[2])
    TD3Pend <- switch(as.character(Scenario), 
                      '0' = LRClipped_L1[2], 
                      '1' = width(GSeq),
                      '2' = LRClipped_L1[2],
                      '3' = LRClipped_G[1])
    
    # Determine which junction is covered by a read
    JunctionCovered <- switch(as.character(Scenario), 
                        '0' = 3, 
                        '1' = 5,
                        '2' = 5,
                        '3' = 3)
    # Determine whether read reaches into 5' or 3' junction
    ClipWithL1 <- c(1:2)[c(blnL1OnLeft, blnL1OnRight)]
    Reaches5P  <- JunctionCovered == 5 | LRClipped_G[ClipWithL1] > MinFullL1
    Reaches3P  <- JunctionCovered == 3 | LRClipped_G[ClipWithL1] > MinFullL1
 
    # Collect info in a vector
    c(LeftClipped_G = LRClipped_G[1], RightClipped_G = LRClipped_G[2], 
      LeftClipped_L1 = LRClipped_L1[1], RightClipped_L1 = LRClipped_L1[2],
      InsPos = InsPos, Scenario = Scenario, TD5Pstart = TD5Pstart, 
      TD5Pend = TD5Pend, TD5Pwidth = TD5Pend - TD5Pstart,
      TD3Pstart = TD3Pstart, TD3Pend = TD3Pend, TD3Pwidth = TD3Pend - TD3Pstart,
      JunctionCovered = JunctionCovered, L1Start = L1Start, L1End = L1End,
      Reaches5P = Reaches5P, Reaches3P = Reaches3P)
  }))
  L1Info <- as.data.frame(L1Info)
  L1Info$Name     <- RLL$GenomeReadList$qname
  L1Info$L1Strand <- 1*blnSameStrand # 0 corresponds to negative strand
  L1Info$L1pos    <- RLL$L1ReadList$pos
  L1Info$Scenario[is.na(L1Info$InsPos)] <- NA
  L1Info
})

# Loop through potential full-length insertions and extract info
for (i in idxPotentialFull){
  
  # Get Info
  L1Info   <- ReadInfoList[[which(i == idxPotentialFull)]]  

  # Fill in info about L1 insertion in FullL1Info
  FullL1Info$L1StrandNum[i] <- mean(L1Info$L1Strand, na.rm = T)
  if (mean(L1Info$L1Strand, na.rm = T) > 0.5) {
    FullL1Info$L1Strand[i] <- "+"
  } else {
    FullL1Info$L1Strand[i] <- "-"
  }
  FullL1Info$L1InsertionPosition.median[i] <- median(L1Info$InsPos)
  FullL1Info$L1InsertionPosition.min[i]    <- min(L1Info$InsPos)
  FullL1Info$L1InsertionPosition.max[i]    <- max(L1Info$InsPos)
  
  # Determine which reads cover which junction
  bln5Covered <- L1Info$JunctionCovered == 5
  bln3Covered <- L1Info$JunctionCovered == 3
  
  # Get start and end of L1 
  FullL1Info$L1Start.median[i] <- median(L1Info$L1Start[L1Info$Reaches5P])
  FullL1Info$L1Start.min[i]    <- min(L1Info$L1Start[L1Info$Reaches5P])
  FullL1Info$L1Start.max[i]    <- max(L1Info$L1Start[L1Info$Reaches5P])
  FullL1Info$L1End.median[i]   <- median(L1Info$L1End[L1Info$Reaches3P])
  FullL1Info$L1End.min[i]      <- min(L1Info$L1End[L1Info$Reaches3P])
  FullL1Info$L1End.max[i]      <- max(L1Info$L1End[L1Info$Reaches3P])

  # Get length of 5' transduced sequence 
  FullL1Info$L15PTransdSeq.median[i] <- median(L1Info$TD5Pwidth[bln5Covered])
  FullL1Info$L15PTransdSeq.min[i]    <- min(L1Info$TD5Pwidth[bln5Covered])
  FullL1Info$L15PTransdSeq.max[i]    <- max(L1Info$TD5Pwidth[bln5Covered])
 
  # Get length of 3' transduced sequence 
  FullL1Info$L13PTransdSeq.median[i] <- median(L1Info$TD3Pwidth[bln3Covered])
  FullL1Info$L13PTransdSeq.min[i]    <- min(L1Info$TD3Pwidth[bln3Covered])
  FullL1Info$L13PTransdSeq.max[i]    <- max(L1Info$TD3Pwidth[bln3Covered])
  
  # Get number of supporting reads
  FullL1Info$NrReads5P[i] <- sum(L1Info$JunctionCovered == 5)
  FullL1Info$NrReads3P[i] <- sum(L1Info$JunctionCovered == 3)
  FullL1Info$NrReadsCover5P[i] <- sum(L1Info$JunctionCovered == 5)
  FullL1Info$NrReadsCover3P[i] <- sum(L1Info$JunctionCovered == 3)
  FullL1Info$NrReadsReach5P[i] <- sum(L1Info$Reaches5P)
  FullL1Info$NrReadsReach3P[i] <- sum(L1Info$Reaches3P)
  FullL1Info$NrSupportReads[i] <- nrow(L1Info)
}

# Write out table with L1 info
write.csv(FullL1Info, file = "D:/L1polymORF/Data/L1InsertionInfo.csv")

# Assess L1 insertion quality
mean(FullL1Info$NrSupportReads)
mean(FullL1Info$NrReadsCover5P)
mean(FullL1Info$NrReadsCover3P)
mean(FullL1Info$NrReadsReach3P)
blnMultipleReads    <- FullL1Info$NrSupportReads > 1
blnConsistentStrand <- (FullL1Info$L1StrandNum == 0 | FullL1Info$L1StrandNum == 1)
blnBothJunctionsCovered <- FullL1Info$NrReads5P > 0 & FullL1Info$NrReads3P > 0
sum(blnMultipleReads & blnConsistentStrand)
sum(blnMultipleReads & blnConsistentStrand & blnBothJunctionsCovered)
sum(blnMultipleReads & blnBothJunctionsCovered)
hist(FullL1Info$NrReads5P)
hist(FullL1Info$NrReads3P)

# Range in estimate of 5' transduced sequence
Range5PTdSeq <- (FullL1Info$L15PTransdSeq.max - FullL1Info$L15PTransdSeq.min)[FullL1Info$NrReadsCover5P > 1]
PropVar5PTdSeq <- Range5PTdSeq / (FullL1Info$L15PTransdSeq.med[FullL1Info$NrReadsCover5P > 1])
hist(PropVar5PTdSeq)
hist(Range5PTdSeq)
hist(FullL1Info$L15PTransdSeq.med)
sum(abs(PropVar5PTdSeq) < 0.1)

###############################################
#                                             #
#       Create alignments of reads            #
#         covering full-length L1             #
#                                             #
###############################################

# Get indices of insertion with strong evidence of full-length insertion and 
# subset info file
rownames(FullL1Info)[which(FullL1Info$Max3P >= 6020)]
idxFull2 <- which(rownames(FullL1Info) %in% FilesWithReads[idxFull])
FullL1InfoSubset <- FullL1Info[idxFull2,]

###############
# Align transduced sequences
###############

# Align transduced sequences and keep info about each insertion
AlignmentFilesTD5P <- c()
AlignmentFilesTD3P <- c()
FastaFiles_5P <- c()
FastaFiles_3P <- c()
ReadInfoListSubset        <- ReadInfoList[idxPotentialFull %in% idxFull2]
names(ReadInfoListSubset) <- FilesWithReads[idxFull]
for (i in idxFull2){
  
  # Get position of current entry in ReadListPerPeak
  x <- idx5P[i]
  
  # Get reads mapped to L1, retain only primary reads and get the number of bp
  # clipped on the left
  RL <- MatchedReadList[[which(i == idxPotentialFull)]]
  
  # Get Info and match reads to info
  L1Info <- ReadInfoList[[which(i == idxPotentialFull)]]  

  # Create a file name for output fasta file
  FastaFile5P <- paste(FastaFilePath, FullL1Info$chromosome[i], "_",
                     FullL1Info$idx[i], "_TDS_5P.fas", sep = "")
  FastaFile3P <- paste(FastaFilePath, FullL1Info$chromosome[i], "_",
                       FullL1Info$idx[i], "_TDS_3P.fas", sep = "")
  FastaFiles_5P <- c(FastaFiles_5P, FastaFile5P)
  FastaFiles_3P <- c(FastaFiles_3P, FastaFile3P)
  
  # Determine which sequence has enough clipped
  idxEnoughClipped <- which(
    pmax(L1Info$LeftClipped_G, L1Info$RightClipped_G) > MinClipAlign &
      L1Info$L1pos < 500)

  # Loop through reads and collect parts of reads that are transduced sequences
  TdS5PList   <- vector(mode = "list", length = length(idxEnoughClipped))
  TdS3PList   <- vector(mode = "list", length = length(idxEnoughClipped))
  idxL5P <- 0
  idxL3P <- 0
  for (y in idxEnoughClipped) {
    
    # Get sequence of current read and the number of left and right clipped bps
    Seq   <- RL$GenomeReadList$seq[y]
    if (L1Info$TD5Pstart[y] < L1Info$TD5Pend[y]){
      TdS5P   <- subseq(Seq, L1Info$TD5Pstart[y], L1Info$TD5Pend[y])
      idxL5P  <- idxL5P + 1
      TdS5PList[[idxL5P]] <- s2c(as.character(TdS5P))
    }
    if (L1Info$TD3Pstart[y] < L1Info$TD3Pend[y]){
      TdS3P   <- subseq(Seq, L1Info$TD3Pstart[y], L1Info$TD3Pend[y])
      idxL3P  <- idxL3P + 1
      TdS3PList[[idxL3P]] <- s2c(as.character(TdS3P))
    }
  }
  
  # Save transduced sequences as fasta file
  names(TdS5PList) <- RL$GenomeReadList$qname[idxEnoughClipped]
  names(TdS5PList) <- gsub("/", "_", names(TdS5PList))
  write.fasta(TdS5PList, names = names(TdS5PList), file.out = FastaFile5P)
  
  names(TdS3PList) <- RL$GenomeReadList$qname[idxEnoughClipped]
  names(TdS3PList) <- gsub("/", "_", names(TdS3PList))
  write.fasta(TdS3PList, names = names(TdS3PList), file.out = FastaFile3P)

  # Aligning transduced sequences as fasta file
  if (length(TdS5PList) > 1){
    AlignedFile5P <- gsub(".fas", "_aligned.fas", FastaFile5P)
    cat("Aligning", FastaFile5P, "\n")
    run_MUSCLE(InputPath = FastaFile5P, OutputPath = AlignedFile5P) 
    AlignmentFilesTD5P <- c(AlignmentFilesTD5P, AlignedFile5P)
  }
  if (length(TdS3PList) > 1){
    AlignedFile3P <- gsub(".fas", "_aligned.fas", FastaFile3P)
    cat("Aligning", FastaFile3P, "\n")
    run_MUSCLE(InputPath = FastaFile3P, OutputPath = AlignedFile3P) 
    AlignmentFilesTD3P <- c(AlignmentFilesTD3P, AlignedFile3P)
  }
}

# Loop through fasta files and create consensus sequences for transduced sequences
TDConsensus_5P <- lapply(FastaFiles_5P, function(x){
  AlignedFile <- gsub(".fas", "_aligned.fas", x)
  if(file.exists(AlignedFile)){
    Alignment <- read.dna(AlignedFile, format = "fasta", as.character = T, as.matrix = T)
    ConsensusSeq <- apply(Alignment, 2, FUN = function(z) {
      Nucs <- z[z != "-"]
      NucFreq <- table(Nucs)
      names(NucFreq)[which.max(NucFreq)]
    })
  } else {
    ConsensusSeq <- read.dna(x, format = "fasta", as.character = T, as.matrix = T)
  }
  ConsensusSeq
})
names(TDConsensus_5P) <- FilesWithReads[idxFull]
TDConsensus_3P <- lapply(FastaFiles_3P, function(x){
  AlignedFile <- gsub(".fas", "_aligned.fas", x)
  if(file.exists(AlignedFile)){
    Alignment <- read.dna(AlignedFile, format = "fasta", as.character = T, as.matrix = T)
    ConsensusSeq <- apply(Alignment, 2, FUN = function(z) {
      Nucs <- z[z != "-"]
      NucFreq <- table(Nucs)
      names(NucFreq)[which.max(NucFreq)]
    })
  } else {
    ConsensusSeq <- read.dna(x, format = "fasta", as.character = T, as.matrix = T)
  }
  ConsensusSeq
})
names(TDConsensus_3P) <- FilesWithReads[idxFull]

# Loop through potential full-length insertions and 
AlignmentFiles <- c()
i <- 68
for (i in idxFull2){
  print(i)
  # Index for transduced sequences
  idxTD <-  which(i == idxFull2)
  
  # Get Info and matched read list
  L1Info <- ReadInfoList[[which(i == idxPotentialFull)]]  
  RL     <- MatchedReadList[[which(i == idxPotentialFull)]]$GenomeReadList 

    # Create a file name for output fasta file
  FastaFile <- paste(FastaFilePath, FullL1Info$chromosome[i], "_",FullL1Info$idx[i],
                     ".fas", sep = "")
  
  # Get insertion location and create a reference sequence with sequence
  # flanking L1 and the consensus L1
  InsLoc  <- FullL1Info$L1InsertionPosition.median[i]
  L1Str   <- FullL1Info$L1Strand[i]
  L1start <- min(L1Info$L1pos)
  L1ConsSeqLocal <- L1CharV[L1start:length(L1CharV)]
  if (L1Str == "-") {
    L1ConsSeqLocal <- L1CharV_RV[1:(length(L1CharV_RV) - L1start + 1)]
    LeftTDs  <-  TDConsensus_3P[[idxTD]]
    RightTDs <-  TDConsensus_5P[[idxTD]]
  } else {
    LeftTDs  <-  TDConsensus_5P[[idxTD]]
    RightTDs <-  TDConsensus_3P[[idxTD]]
  }
  TransDL    <- FullL1Info$L15PTransdSeq.median[i]
  LeftFlankR <- GRanges(FullL1Info$chromosome[i], 
                        IRanges(start = InsLoc - FlankSize + 1,
                                   end = InsLoc))
  RightFlankR <- GRanges(FullL1Info$chromosome[i], IRanges(start = InsLoc + 1,
                            end = InsLoc + FlankSize))
  LeftFlankSeq    <- getSeq(BSgenome.Hsapiens.UCSC.hg19, LeftFlankR)
  RightFlankSeq   <- getSeq(BSgenome.Hsapiens.UCSC.hg19, RightFlankR)
  LeftFlankCharv  <- s2c(as.character(LeftFlankSeq))
  RightFlankCharv <- s2c(as.character(RightFlankSeq))
  
  RefSeqCharV     <- c(LeftFlankCharv, "N", LeftTDs, "N", 
                       L1ConsSeqLocal, "N", RightTDs, "N", RightFlankCharv)
  RefSeqCharV     <- toupper(RefSeqCharV)
                       
  # Get reads mapped to the genome on current locus
  LRClippedG <- sapply(RL$cigar, NrClippedFromCigar)
  idxEnoughClipped <- which(
    pmax(LRClippedG[1,], LRClippedG[2,]) > MinClipAlign)

  # Loop through reads and collect parts of reads that should align with L1 and
  # its flanking sequence
  SeqList   <- vector(mode = "list", length = length(idxEnoughClipped) + 1)
  SeqList[[1]]     <- RefSeqCharV
  idxL <- 1
  for (y in idxEnoughClipped) {
    
    # Get sequence of current read and the number of left and right clipped bps
    Seq        <- RL$seq[y]
    LRClipped  <- LRClippedG[,y]
    
    # If more is clipped on the left take the sequence from the start into the
    # flank extending to the right of the cilpper part
    SeqStart <- switch(as.character(L1Info$Scenario[y]), 
                        '0' = max(1, L1Info$TD3Pstart[y] - FlankSize), 
                        '1' = max(1, L1Info$TD5Pstart[y] - FlankSize),
                        '2' = 1,
                        '3' = 1)
    SeqEnd <- switch(as.character(L1Info$Scenario[y]), 
                       '0' = width(Seq), 
                       '1' = width(Seq),
                       '2' = min(width(Seq), L1Info$TD5Pend[y] + FlankSize),
                       '3' = min(width(Seq), L1Info$TD3Pend[y] + FlankSize))
    Seq2Add  <- subseq(Seq[[1]], SeqStart, SeqEnd)
    idxL <- idxL + 1
    SeqList[[idxL]] <- s2c(as.character(Seq2Add))
  }
  names(SeqList) <- c("RefSeq", RL$qname[idxEnoughClipped])
  names(SeqList) <- gsub("/", "_", names(SeqList))
  write.fasta(SeqList, names = names(SeqList), file.out = FastaFile)
  
  AlignedFile <- gsub(".fas", "_aligned.fas", FastaFile)
  cat("Aligning", FastaFile, "\n")
  run_MUSCLE(InputPath = FastaFile, OutputPath = AlignedFile) 
  
  # Read alignment and reorder so that reference sequence is first
  Alignment <- read.fasta(AlignedFile)
  idxRefSeq <- which(names(Alignment) == "RefSeq")
  AlignmentReordered <- Alignment
  AlignmentReordered[[1]] <- Alignment[[idxRefSeq]]
  names(AlignmentReordered)[1] <- names(Alignment)[idxRefSeq]
  AlignmentReordered[[idxRefSeq]] <- Alignment[[1]]
  names(AlignmentReordered)[idxRefSeq] <- names(Alignment)[1]
  write.fasta(AlignmentReordered, names = names(AlignmentReordered), file.out = AlignedFile)
  
  # Keep track of alignment file names
  AlignmentFiles <- c(AlignmentFiles, AlignedFile)
}

# Loop through alignment files and analyze regions of 
AlignFile <- AlignmentFiles[1]
AlignInfoList <- lapply(AlignmentFiles, function(AlignFile){
  
  # Extract plot title from file name
  Fsplit1 <- strsplit(AlignFile, "/")[[1]]
  Fsplit2 <- strsplit(Fsplit1[length(Fsplit1)], "_")[[1]]
  Title <- paste(Fsplit2[1:2], collapse = "-")
  
  # Read alignment and extract infor about alignment quality
  Alignment <- read.dna(AlignFile, format = "fasta", as.character = T, as.matrix = T)
  blnSame   <- sapply(1:ncol(Alignment), function(x){
    any(Alignment[2:nrow(Alignment),x] == Alignment[1,x])
  })
  blnSameTS <- as.ts(1*blnSame)
  SmoothedSame <- sma(blnSameTS, order = 20)
  nPos <- which(Alignment[1,] %in% c("n", "N"))
  
  # Extract consensus sequence
  Alignment <- Alignment[ , Alignment[1,] != "-"]
  # ConsensusSeq <- apply(Alignment, 2, FUN = function(z) {
  #   Nucs <- z[z != "-"]
  #   if (z[1] != "n"){
  #     NucFreq <- table(Nucs)
  #     names(NucFreq)[which.max(NucFreq)]
  #   } else {
  #     z[1]
  #   }
  # })
  ConsensusSeq <- sapply(1:ncol(Alignment), function(v) {
    z <- Alignment[,v]
    Nucs <- z[z != "-"]
    if (z[1] != "n" & v > min(nPos) & v < max(nPos)){
      NucFreq <- table(Nucs)
      names(NucFreq)[which.max(NucFreq)]
    } else {
      z[1]
    }
  })
  ConsensNPos <- which(ConsensusSeq %in% c("n", "N"))
  
  
  list(nPos = nPos, SmoothedSame = SmoothedSame, Title = Title,
       ConsensusSeq = ConsensusSeq, ConsensNPos = ConsensNPos)
})

# Get consensus sequence of flanking sequence and insertion
L1WithFlank <- lapply(AlignInfoList, function(x) x$ConsensusSeq)
names(L1WithFlank) <- sapply(AlignInfoList, function(x) x$Title)
names(L1WithFlank) <- gsub("-", "_", names(L1WithFlank))
write.fasta(L1WithFlank, names = names(L1WithFlank),
            file.out = L1InsertionFastaPath)

# Write reconstructed L1 insertions as individual fasta files so that reads
# can be aligned to them
i <- 4
for (i in 1:length(L1WithFlank)){
  FileNameL1 <- paste(FastaFilePath, "L1WithFlank", FlankSize, 
                  "bp_",names(L1WithFlank)[i], ".fas", sep = "")
  FileNameTdSeq <- paste(FastaFilePath, "TdSeq", FlankSize, 
                      "bp_",names(L1WithFlank)[i], ".bed", sep = "")
  FileNameL1pos <- paste(FastaFilePath, "L1pos", FlankSize, 
                         "bp_",names(L1WithFlank)[i], ".bed", sep = "")
  write.fasta(L1WithFlank[[i]], names = names(L1WithFlank)[i], FileNameL1)
  NPos      <- AlignInfoList[[i]]$ConsensNPos
  TransdPos <- data.frame(rep(names(L1WithFlank)[i], 2), NPos[c(1,3)], 
                          NPos[c(2,4)])
  L1Pos <- data.frame(names(L1WithFlank)[i], NPos[2], NPos[3])
  write.table(TransdPos, file = FileNameTdSeq, col.names = F, quote = F, 
              row.names = F, sep = "\t")
  write.table(L1Pos, file = FileNameL1pos, col.names = F, quote = F, 
              row.names = F, sep = "\t")
}

# Plot alignment quality
par(mfrow = c(3, 2), mar = c(3, 2, 3, 0.5), oma = c(2, 4, 1, 1))
for (i in 1:length(AlignInfoList)){
  SmoothedSame <- AlignInfoList[[i]]$SmoothedSame
  nPos <- AlignInfoList[[i]]$nPos
  Title <- AlignInfoList[[i]]$Title
  plot(SmoothedSame$fitted, type = "l", xlab = "", ylab = "", main = Title,
       ylim = c(0, 1.1))
  segments(x0 = nPos, y0 = 0, y1 = 1, col = "red", lwd = 2)
}
mtext("Proportion of bp same as in reference", 2, outer = T, line = 2)
CreateDisplayPdf("D:/L1polymORF/Figures/InsertionAlignment.pdf",
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 10, width = 7)

# Get list of sequences with insertions
InsertionSeqs <- lapply(idxFull2, function(i){
  
  # Get Info
  L1Info <- ReadInfoList[[which(i == idxPotentialFull)]]  
  
  # Get insertion location and create a reference sequence with sequence
  # flanking L1 and the consensus L1
  InsLoc     <- FullL1Info$L1InsertionPosition.median[i]
  L1Str      <- NewStrand[i == idxFull2]
  L1ConsSeqLocal <- L1CharV
  if (L1Str == "-") L1ConsSeqLocal <- L1CharV_RV
  LeftFlankR <- GRanges(FullL1Info$chromosome[i], 
                        IRanges(start = InsLoc - FlankSize + 1,
                                end = InsLoc))
  RightFlankR <- GRanges(FullL1Info$chromosome[i], IRanges(start = InsLoc + 1,
                                                           end = InsLoc + FlankSize))
  LeftFlankSeq    <- getSeq(BSgenome.Hsapiens.UCSC.hg19, LeftFlankR)
  RightFlankSeq   <- getSeq(BSgenome.Hsapiens.UCSC.hg19, RightFlankR)
  LeftFlankCharv  <- s2c(as.character(LeftFlankSeq))
  RightFlankCharv <- s2c(as.character(RightFlankSeq))
  c(LeftFlankCharv, L1ConsSeqLocal, RightFlankCharv)
})
SeqNames <- paste(FullL1InfoSubset$chromosome, FullL1InfoSubset$idx, sep = "_")
write.fasta(InsertionSeqs, names = SeqNames, file.out = "D:/L1polymORF/Data/L1InsertionSequences.fas")
InsertionSummaries <- FullL1InfoSubset[,c("chromosome", "L1InsertionPosition.median")]
colnames(InsertionSummaries)[colnames(InsertionSummaries) == "L1InsertionPosition.median"] <- "L1InsertionPosition"
InsertionSummaries$SeqName  <- SeqNames
InsertionSummaries$L1Strand <- NewStrand
InsertionSummaries$TransD5P <- c("Yes", "No", "No", "Yes",  "Yes")
InsertionSummaries$TransD3P <- c("No", "Yes", "No", "No",  "No")
write.csv(InsertionSummaries, file = "D:/L1polymORF/Data/L1InsertionSummary.csv",
          row.names = F)

# Create genomic ranges for L1 insertions and compare with known L1 ranges
GRL1Capture <- makeGRangesFromDataFrame(FullL1Info)
GRL1Capture100 <- resize(GRL1Capture, 100, fix = "center")
L1CatalogGR100 <- resize(L1CatalogGR, 100, fix = "center")
GRL1Ins1000G100 <- resize(GRL1Ins1000G, 100, fix = "center")
if (any(c(sum(overlapsAny(GRL1Capture100, L1CatalogGR100)),
          sum(overlapsAny(GRL1Capture100, GRL1Ins1000G100))) > 0)){
  cat("Some newly found L1 insertion overlap with catalog or 1000 Genome elements!\n")
}

###############################################
#                                             #
#       Assemble sequences of L1              #
#             insertions                      #
#                                             #
###############################################

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

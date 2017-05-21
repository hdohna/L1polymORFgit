# The script below reads a list of reads aligned to L1, a coverage and quantile 
# matrix (created in script 'CalcCoverMatReadList.R') and creates new L1 
# reference with known and new L1 stumps

##############
# Source prerequisites
##############

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(seqinr)
library(ShortRead)
library(Rsamtools)
library(rtracklayer)
library(csaw)
library(GenomicRanges)
#library(BSgenome.Hsapiens.UCSC.hg19)

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
MaxClip    <- 100
MinL1End   <- 6000
MaxL1Start <- 100

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
# Path to L1 GRanges from 1000 genome data
CoveragePlotPath      <- 'D:/L1polymORF/Figures/L1InsertionCoverage_NA12878_PacBio.pdf'
CoverDataPath         <- 'D:/L1polymORF/Data/L1_NA12878_PacBio_Coverage.RData'
CoverDataPath         <- 'D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_Results.Rdata'
CoverDataPath         <- 'D:/L1polymORF/Data/PacBioHiFi__Results.Rdata'
PeakRangePath         <- "D:/L1polymORF/Data/PacBioHiFiSubreads_L1Ranges.RData"
RefRangePath          <- "D:/L1polymORF/Data/L1RefRanges_hg19.RData"
L1_1000GenomeDataPath <- "D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData"
L1_1000_NA12878_Path  <- "D:/L1polymORF/Data/NA12878.1000genome.L1HS.insert.bed"
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
load(RefRangePath)
load(PeakRangePath)
load(L1_1000GenomeDataPath)

# Load coverage data (calculated in CalcCoverMatReadList)
load(CoverDataPath)

######
# Read and process L1 catalog
######

# Read in table with known L1 
L1Catalogue <- read.csv("D:/L1polymORF/Data/L1CatalogExtended.csv", as.is = T)

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
L1CatalogL1Mapped <- L1Catalogue[which(blnL1Mapped & blnAllele1),]

# Lift coordinates and get genomic ranges for catalog L1 on hg19
LiftOverList <- LiftoverL1Catalog(L1CatalogL1Mapped,
                                  ChainFilePath = "D:/L1polymORF/Data/hg38ToHg19.over.chain")
L1CatalogGR <- LiftOverList$GRCatalogue_hg19# Specify folder 

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
RepeatTable <- read.csv("D:/L1polymORF/Data/repeatsHg19_L1HS.csv")
#RepeatTable <- RepeatTable[RepeatTable$repName == "L1HS", ]

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
L1RefGR <- GRanges(seqnames = RepeatTable$genoName,
                   ranges = IRanges(start = RepeatTable$genoStart,
                                    end = RepeatTable$genoEnd),
                   strand = RepeatTable$strand)
L1RefGRFull <- L1RefGR[width(L1RefGR) > 6000]

# L1RefSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, L1RefGRFull)
# L1RefSeq <- as.character(L1RefSeq)
# StrandV  <- as.vector(strand(L1RefGRFull))
# L1RefSeqMat <- sapply(1:length(L1RefSeq), function(i) {
#   L1V <- strsplit(L1RefSeq[i], "")[[1]]
#   if (StrandV[i] == "+"){
#     L1V[1:FivePEnd]
#   } else {
#     L1V[(length(L1V) - FivePEnd + 1):length(L1V)]
#   }
# })
# dim(L1RefSeqMat)
# colnames(L1RefSeqMat) <- paste(as.vector(seqnames(L1RefGRFull)), start(L1RefGRFull),
#                                end(L1RefGRFull), "Ref", sep = "_")

###############################################
#                                             #
#  Collect info on promising L1 insertions    #
#            (FullL1Info)                     #
#                                             #
###############################################

# Collect information on insertion that fullfill a certain minimum criterion
idx5P <- which(sapply(1:nrow(CoverMat), function(x) {
  all(CoverMat[x, FivePStart:FivePEnd] >= MinCover)
}))
cat("******  ", length(idx5P), "insertions have a minimum coverage of", 
MinCover, "in 5' region from", FivePStart, "to", FivePEnd, "   *********\n")
idxFull <- which(sapply(1:nrow(CoverMat), function(x) all(CoverMat[x, FivePStart:6040] > 0)))
FullL1Info <- t(sapply(FilesWithReads, function(x){
  FPathSplit <- strsplit(x, "/")[[1]]
  FName      <- FPathSplit[length(FPathSplit)]
  FName      <- substr(FName, 1, nchar(FName) - 4)
  strsplit(FName, "_")[[1]]
}))
FullL1Info <- base::as.data.frame(FullL1Info, stringsAsFactors = F)
colnames(FullL1Info) <- c("chromosome", "idx")
FullL1Info$Max3P <- sapply(1:nrow(CoverMat), function(x) max(which(CoverMat[x, ] > 0)))
FullL1Info$idx   <- as.numeric(FullL1Info$idx)
FullL1Info$start <- start(IslGRanges_reduced[FullL1Info$idx])
FullL1Info$end   <- end(IslGRanges_reduced[FullL1Info$idx])
FullL1Info$cover   <- CoverMat[, 100]
FullL1Info$First0Cover <- sapply(1:nrow(CoverMat), 
   function(x) min(50 + which(CoverMat[x, 50:ncol(CoverMat)] == 0)))
FullL1Info$FirstCover <- sapply(1:nrow(CoverMat), 
   function(x) min(which(CoverMat[x, ] > 0)))

idxFull2 <- which(rownames(FullL1Info) %in% FilesWithReads[idxFull])

# Set up insertion descriptors
FullL1GR <- makeGRangesFromDataFrame(FullL1Info)
FullL1Info$PotentialFullLength        <- NA
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
FullL1Info$NrZMW_ID                   <- NA
FullL1Info$Dist2_1000G                <- NA
FullL1Info$Dist2_NA12878              <- NA
FullL1Info$Dist2_Catalog              <- NA
FullL1Info$L1LengthClosest_1000G      <- NA
FullL1Info$L1LengthClosest_NA12878    <- NA

# Loop through list of matching reads and collect info
ReadInfoList <- lapply(1:length(MatchedReadList), function(x) {
  print(x)
  RLL <- MatchedReadList[[x]]
  # Determine whether reads are on the same strand  
  blnSameStrand <- RLL$GenomeReadList$strand == RLL$L1ReadList$strand
  
  # Loop through reads an collect information on number of bps clipped on
  # reads mapped to genome, the insertion location and the insertion strand
  ClipWithL1 <- 0
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
    MotifStart <- min(LRClipped_L1[1] + MotifOffSet, width(L1Seq))
    MotifEnd   <- min(LRClipped_L1[1] + MotifOffSet + MotifWidth, width(L1Seq))
    L1Motif    <- subseq(L1Seq, MotifStart, MotifEnd)
    if (width(L1Motif) < MotifWidth){
      warning("Motif width is reduced to ", width(L1Motif), "\n")
    }
    if (!blnSameStrand[y]) {L1Motif <- reverseComplement(L1Motif)}
    motifMatchLeft  <- matchPattern(L1Motif[[1]], LeftClippedSeq[[1]])
    motifMatchRight <- matchPattern(L1Motif[[1]], RightClippedSeq[[1]])
    blnL1OnLeft     <- length(motifMatchLeft)  > 0
    blnL1OnRight    <- length(motifMatchRight) > 0
    
    # Identify position of L1 insertion
    InsPos <- NA
    if (blnL1OnLeft){
      InsPos <- RLL$GenomeReadList$pos[y]
    } 
    if (blnL1OnRight){
      InsPos <- RLL$GenomeReadList$pos[y] + 
        ReadLengthFromCigar(RLL$GenomeReadList$cigar[y])
    } 
    # Scenario: 0 = L1 on negative strand and read maps on left side of L1
    #           1 = L1 on positive strand and read maps on left side of L1
    #           2 = L1 on negative strand and read maps on right side of L1
    #           3 = L1 on positive strand and read maps on right side of L1
    Scenario <- blnSameStrand[y] + 2 * blnL1OnLeft
    
    if (blnL1OnLeft & blnL1OnRight){
      warning("L1 motif matches left and right clipped sequence in", y, "\n", immediate. = T)
      Scenario   <- NA
      ClipWithL1 <- NA
    }
    if ((length(motifMatchLeft) == 0) & (length(motifMatchRight) == 0)){
      #browser()
      warning("L1 motif matches no clipped sequence in read", y, "\n")
      Scenario   <- NA
      ClipWithL1 <- NA
    }

    # Get start and stop indices of transduced sequence of 5' and 3' end
    TD5Pstart <- switch(as.character(Scenario), 
                        '0' = width(GSeq) - LRClipped_L1[1], 
                        '1' = width(GSeq) - LRClipped_G[2],
                        '2' = width(GSeq) - LRClipped_L1[1],
                        '3' = 1, NA)
    TD5Pend   <- switch(as.character(Scenario), 
                        '0' = width(GSeq), 
                        '1' = LRClipped_L1[1],
                        '2' = LRClipped_G[1],
                        '3' = LRClipped_L1[1], NA)
    TD3Pstart <- switch(as.character(Scenario), 
                        '0' = width(GSeq) - LRClipped_G[2], 
                        '1' = width(GSeq) - LRClipped_L1[2],
                        '2' = 1,
                        '3' = width(GSeq) - LRClipped_L1[2], NA)
    TD3Pend <- switch(as.character(Scenario), 
                      '0' = LRClipped_L1[2], 
                      '1' = width(GSeq),
                      '2' = LRClipped_L1[2],
                      '3' = LRClipped_G[1], NA)
    
    # Determine which junction is covered by a read
    JunctionCovered <- switch(as.character(Scenario), 
                        '0' = 3, 
                        '1' = 5,
                        '2' = 5,
                        '3' = 3, NA)
    # Determine whether read reaches into 5' or 3' junction
    if (!is.na(ClipWithL1)){
      ClipWithL1 <- c(1:2)[c(blnL1OnLeft, blnL1OnRight)]
      Reaches5P  <- JunctionCovered == 5 | LRClipped_G[ClipWithL1] > MinFullL1
      Reaches3P  <- JunctionCovered == 3 | LRClipped_G[ClipWithL1] > MinFullL1
    } else {
      Reaches5P  <- NA
      Reaches3P  <- NA
    }
    
    # Get ZMW id from read ID
    ZMW_ID <- strsplit(RLL$L1ReadList$qname[y], "/")[[1]][2]
    ZMW_ID <- as.numeric(ZMW_ID)
    
    # Collect info in a vector
    InfoV <- c(LeftClipped_G = LRClipped_G[1], RightClipped_G = LRClipped_G[2], 
      LeftClipped_L1 = LRClipped_L1[1], RightClipped_L1 = LRClipped_L1[2],
      InsPos = InsPos, Scenario = Scenario, TD5Pstart = TD5Pstart, 
      TD5Pend = TD5Pend, TD5Pwidth = TD5Pend - TD5Pstart,
      TD3Pstart = TD3Pstart, TD3Pend = TD3Pend, TD3Pwidth = TD3Pend - TD3Pstart,
      JunctionCovered = JunctionCovered, L1Start = L1Start, L1End = L1End,
      Reaches5P = Reaches5P, Reaches3P = Reaches3P, ZMW_ID = ZMW_ID)
    if (length(InfoV) != 18) browser()
    InfoV
  }))
  L1Info          <- as.data.frame(L1Info)
  L1Info$Name     <- RLL$GenomeReadList$qname
  L1Info$L1Strand <- 1*blnSameStrand # 0 corresponds to negative strand
  L1Info$L1pos    <- RLL$L1ReadList$pos
  L1Info$Scenario[is.na(L1Info$InsPos)] <- NA
  L1Info
})

# Loop through insertions and extract info
for (i in 1:length(ReadInfoList)){
  
  # Get Info
  L1Info   <- ReadInfoList[[i]]  
  
  # Check whether any read has a long left-clipped pr right-clipped sequence 
  # on L1 that indicates L1 insertion might not be full-length
  blnNotFull <- (L1Info$LeftClipped_L1 > MaxClip & L1Info$L1Start > MaxL1Start) |
                (L1Info$RightClipped_L1 > MaxClip & L1Info$L1End < MinL1End)
  FullL1Info$PotentialFullLength[i] <- !blnNotFull
  
  # Fill in info about L1 insertion (strandedness and position)
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
  FullL1Info$NrZMW_ID[i]       <- length(unique(L1Info$ZMW_ID))
}
cat("******  ", sum(FullL1Info$PotentialFullLength, na.rm = T), 
    "insertions are potentially full length *********\n")
warnings()

# Write out table with L1 info
write.csv(FullL1Info, file = "D:/L1polymORF/Data/L1InsertionInfo.csv")


###############################################
#                                             #
#       Compare to 1000 genome data           #
#                                             #
###############################################

# Create genomic ranges from table 
L1InsGR <- makeGRangesFromDataFrame(FullL1Info)

# Read in 1000 Genome L1 insertion (created in script Create_1000G_L1GRanges.R)
load(L1_1000GenomeDataPath)

# Calculate distance from each L1 to closest L1 in 1000 genome
Dist2_1000GL1 <- distanceToNearest(L1InsGR, L1_1000G_GR_hg19, ignore.strand = T) 
FullL1Info$Dist2_1000G[Dist2_1000GL1@from]  <- Dist2_1000GL1@elementMetadata@listData$distance

# Record the insertion length of closest L1 in 1000 genome
L1Length <- L1_1000G_GR_hg19@elementMetadata@listData$InsLength
FullL1Info$L1LengthClosest_1000G[Dist2_1000GL1@from] <- L1Length[Dist2_1000GL1@to]

# Calculate distance from each L1 to closest NA1287 L1 in 1000 genome
Dist2_1000GL1 <- distanceToNearest(L1InsGR, L1_1000G_GR_hg19_NA12878, ignore.strand = T) 
FullL1Info$Dist2_NA12878[Dist2_1000GL1@from]  <- Dist2_1000GL1@elementMetadata@listData$distance

# Record the insertion length of closest L1 in 1000 genome
L1Length_NA12878 <- L1_1000G_GR_hg19_NA12878@elementMetadata@listData$InsLength
FullL1Info$L1LengthClosest_NA12878[Dist2_1000GL1@from] <- L1Length_NA12878[Dist2_1000GL1@to]

# Calculate distance from each L1 to closest NA1287 L1 in catalog
idxCat_NA12878 <- which(L1CatalogL1Mapped$Coriell_ID == "NA12878")
Dist2_Cat <- distanceToNearest(L1InsGR, L1CatalogGR, ignore.strand = T) 
FullL1Info$Dist2_Catalog[Dist2_Cat@from]  <- Dist2_Cat@elementMetadata@listData$distance

# Look at catalog elements that are colse
idxCatClose <- Dist2_Cat@to[which(Dist2_Cat@elementMetadata@listData$distance < 50)]
L1CatalogL1Mapped[idxCatClose, c("Accession", "Allele", "Chromosome", "Activity", 
   "Allele_frequency", "Reference", "Coriell_ID")]


###############################################
#                                             #
#      Generate a report                      #
#                                             #
###############################################

# Plot average coverage
par(mfrow = c(1,1))
plot(colMeans(CoverMat), type = "l", xlab = "Position on L1", ylab = "Average coverage")
CreateDisplayPdf("D:/L1polymORF/Figures/AverageCoverPerL1.pdf",
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7, width = 7)

# Plot first position where coverage ends
hist(FullL1Info$First0Cover, xlab = "Position on L1 where coverage ends")
hist(FullL1Info$FirstCover, xlab = "Position on L1 where coverage starts")

# Create boolean vectors of various indicators
blnCloseNA12878     <- FullL1Info$Dist2_NA12878 < 50
blnCloseNA12878Full <- FullL1Info$L1LengthClosest_NA12878 > 6000
blnClose1000G       <- FullL1Info$Dist2_1000G < 50
blnClose1000GFull   <- FullL1Info$L1LengthClosest_1000G > 6000
blnAllCovered       <- sapply(1:nrow(CoverMat), function(x) all(CoverMat[x, 100:5000] > 0))
blnMultipleReads    <- FullL1Info$NrSupportReads > 1
blnMultipleZMW      <- FullL1Info$NrZMW_ID > 1
blnConsistentStrand <- (FullL1Info$L1StrandNum == 0 | FullL1Info$L1StrandNum == 1)
blnBothJunctionsCovered <- FullL1Info$NrReads5P > 0 & FullL1Info$NrReads3P > 0

# Function that determines the percentage of full-length insertions for a given
#  filter
Nr1000G <- function(Filter){
  sum(blnClose1000G[Filter], na.rm = T)
}
PropFullL1 <- function(Filter){
  sum(blnClose1000G[Filter] & blnClose1000GFull[Filter], na.rm = T) /
    sum(blnClose1000G[Filter], na.rm = T)
}

# Plot how the proportion of full-length L1 changes with different 
# thresold criteria
StartThreshold <- seq(100, 6000, 100)
PropFullPerStart <- sapply(StartThreshold, function(x) {
  xFilter <- FullL1Info$FirstCover >= x
  PropFullL1(xFilter)
})
plot(StartThreshold, PropFullPerStart, type= "s", 
     xlab = "Coverage on L1 starts after", ylab = "Proportion of full-length L1")
PropFullPerEnd <- sapply(StartThreshold, function(x) {
  xFilter <- FullL1Info$First0Cover >= x
  PropFullL1(xFilter)
})
NrPerEnd <- sapply(StartThreshold, function(x) {
  xFilter <- FullL1Info$First0Cover >= x
  Nr1000G(xFilter)
})
plot(StartThreshold, PropFullPerEnd, type= "s",
     xlab = "Coverage on L1 ends after", ylab = "Proportion of full-length L1")
CreateDisplayPdf("D:/L1polymORF/Figures/PropFullL1.pdf",
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7, width = 7)

# Range in estimate of 5' transduced sequence
Range5PTdSeq <- (FullL1Info$L15PTransdSeq.max - FullL1Info$L15PTransdSeq.min)[FullL1Info$NrReadsCover5P > 1]
PropVar5PTdSeq <- Range5PTdSeq / (FullL1Info$L15PTransdSeq.med[FullL1Info$NrReadsCover5P > 1])
hist(PropVar5PTdSeq)
hist(Range5PTdSeq)
hist(FullL1Info$L15PTransdSeq.med)
sum(abs(PropVar5PTdSeq) < 0.1)

# Create a table with basic quantities
SummaryReport <- data.frame(Description = "Peaks with at least 1 read mapped to L1", 
                              Count = nrow(CoverMat), 
                              stringsAsFactors = F)
SummaryReport <- rbind(SummaryReport,
                         c("L1 insertions supported by > 1 ZMW",
                           sum(blnMultipleZMW & blnConsistentStrand, na.rm = T)))
SummaryReport <- rbind(SummaryReport,
                         c("L1 insertions overlapping with 1000 Genome L1s",
                           sum(blnClose1000G, na.rm = T)))
SummaryReport <- rbind(SummaryReport,
                         c("L1 insertions overlapping with 1000 Genome L1s of NA12878",
                           sum(blnCloseNA12878, na.rm = T)))
SummaryReport <- rbind(SummaryReport,
   c("L1 insertions supported by > 1 ZMW & overlapping with 1000 Genome L1s",
     sum(blnMultipleZMW & blnConsistentStrand & blnClose1000G, na.rm = T)))
SummaryReport <- rbind(SummaryReport,
   c("L1 insertions overlapping with 1000G L1s & L1 > 6000 bp",
      sum(blnClose1000G & blnClose1000GFull, na.rm = T)))
SummaryReport <- rbind(SummaryReport,
    c("L1 insertions overlapping with 1000G L1s of NA12878 & L1 > 6000 bp",
          sum(blnCloseNA12878 & blnCloseNA12878Full, na.rm = T)))
SummaryReport <- rbind(SummaryReport,
   c("Proportion of full-length 1000G L1s of NA12878 detected by PacBio",
     sum(blnCloseNA12878 & blnCloseNA12878Full, na.rm = T) / sum(L1Length_NA12878 > 6000, na.rm = T)))
SummaryReport <- rbind(SummaryReport,
    c("Proportion of L1 with > 1 ZMW among full-length 1000G L1s of NA12878 detected by PacBio",
    sum(blnMultipleZMW & blnCloseNA12878 & blnCloseNA12878Full, na.rm = T)/ sum(blnCloseNA12878 & blnCloseNA12878Full, na.rm = T)))

write.csv(SummaryReport, "D:/L1polymORF/Data/PacBioL1Calls.csv", row.names = F)

# Assess L1 insertion quality
sum(blnAllCovered)
mean(FullL1Info$NrSupportReads)
sum(FullL1Info$NrZMW_ID > 1 & FullL1Info$PotentialFullLength, na.rm = T)
mean(FullL1Info$NrReadsCover5P, na.rm = T)
mean(FullL1Info$NrReadsCover3P, na.rm = T)
mean(FullL1Info$NrReadsReach3P, na.rm = T)
sum(blnMultipleReads)
sum(blnMultipleReads & blnConsistentStrand)
sum(blnMultipleReads & blnConsistentStrand & blnBothJunctionsCovered, na.rm = T)
sum(blnMultipleReads & blnBothJunctionsCovered)
hist(FullL1Info$NrReads5P)
hist(FullL1Info$NrReads3P)
which(idx5P %in% which(blnAllCovered))
rownames(CoverMat)[blnAllCovered]
which(blnMultipleReads & blnBothJunctionsCovered)
FullL1Info[blnMultipleReads & blnBothJunctionsCovered, c("chromosome", "idx")]
FullL1Info$L1InsertionPosition.median[blnMultipleReads & blnBothJunctionsCovered]

hist(FullL1Info$Max3P)
hist(FullL1Info$First0Cover)

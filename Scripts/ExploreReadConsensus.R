# This script explores aspects of read lists

# Load packages
library(seqinr)
library(ShortRead)
library(Rsamtools)
library(rtracklayer)
library(csaw)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Matrix)

# Paths (local computer)
PathStart <- 'D:/L1polymORFgit/Scripts/_Start_L1polymORF.R'
PathBam   <- "D:/L1polymORF/Data/BZ_NonRef/chr16_4873.bam"
PathRef   <- "D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fas"

# Source start script
source(PathStart)

# Load data on chromosme length
load("D:/L1polymORF/Data/ChromLengthsHg19.Rdata")

# Read example Bam file
ReadList <- scanBam(PathBam)

# Read in reference sequence
RefSeqAll <- read.fasta(PathRef)
RefSeqAll <- toupper(RefSeqAll$L1HS_L1_Homo_sapien)

# Create genomic ranges
#GR <- GRanges('L1HS_L1_Homo_sapiens', IRanges(10, 1000))
Chrom <- 'chr1'
ChrL   <- ChromLengthsHg19[Chrom]
ChromStarts <- c(seq(1, ChrL, 10^6), ChrL)
GR <- GRanges('chr1', IRanges(ChromStarts[-length(ChromStarts)], ChromStarts[-1]))

# Remove reads with NA positions
RL <- ReadList[[1]]
blnPosNotNA <- !is.na(RL$pos)
RL          <- lapply(RL, function(x) x[blnPosNotNA])

# Get ZMW id from read ID
ZMW_IDs <- sapply(1:length(RL$pos), function(j) {
  strsplit(RL$qname[j], "/")[[1]][2]
})

# Test correlation between phred and mismatch
RL$seq
ReadGR   <- ReadList2GRanges(RL)
CodeV  <- c()
PhredV <- c()
i <- 1
for (i in 1:100){
  print(i)
  SeqV    <- SeqFromCigar(RL$cigar[i], RL$seq[i])
  Phreds  <- PhredFromCigar(RL$cigar[i], RL$qual[i])
  if (length(Phreds) != length(SeqV)) browser()
  RefSeq  <- getSeq(BSgenome.Hsapiens.UCSC.hg19, ReadGR[i])
  RefSeqV <- strsplit(as.character(RefSeq), "")[[1]]
  RefSeqV <- RefSeqV[-length(RefSeqV)]
  CodeV   <- c(CodeV,  SeqV == RefSeqV)
  PhredV  <- c(PhredV, Phreds)
}
PhredAgg <- aggregate(CodeV, by = list(PhredV), FUN = mean)
plot(PhredAgg)
length(CodeV)
length(PhredV)
# Test various functions

SeqMat   <- SeqMatFromReads(RL, GR[1])
PhredMat <- PhredMatFromReads(RL, GR)
ConsSeq  <- ConsensusFromReads(RL, GR)
RefSeq   <- RefSeqAll[start(GR):end(GR)]
hist(PhredMat[PhredMat > -1])
subseq(RL$seq[1], 59, 61)
# Get a code vector
PredictR <- 2:(nrow(SeqMat) - 1)
CodeV <- c()
FeatureMat <- matrix(ncol = 1, nrow = 0)
for (i in 1:ncol(SeqMat)){
  CodeV <- c(CodeV, SeqMat[PredictR, i] == RefSeq[PredictR])
  FeatureMat <- c(FeatureMat, PhredMat[PredictR, i])
}

# Create a feature matrix
Triplets <- paste(SeqMat[PredictR -1, 1], SeqMat[PredictR, 1], SeqMat[PredictR + 1, 1],
                    sep = "")
UniqueTripl <- unique(Triplets)  
FeatureMat <- sparseMatrix(i = 1:length(Triplets), j = match(Triplets, UniqueTripl))
# FeatureMat <- cbind(FeatureMat, PhredMat[PredictR, 1])
# FeatureMat[1:3, 114]

ME <- maxent(feature_matrix = FeatureMat, code_vector = CodeV, verbose = T)
predict.maxent(ME, feature_matrix = FeatureMat)

PhredVals <- PhredMat[PredictR, 1]
PByPhred <- aggregate(CodeV, by = list(FeatureMat), FUN = mean)
plot(PByPhred)

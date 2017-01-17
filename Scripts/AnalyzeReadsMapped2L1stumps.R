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

# Set file paths
CoveragePlotPath      <- 'D:/L1polymORF/Figures/L1InsertionCoverage_NA12878_PacBio.pdf'
CoverDataPath         <- 'D:/L1polymORF/Data/L1_NA12878_PacBio_Coverage.RData'
L1_1000GenomeDataPath <- "D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData"
OutFolderName_NonRef  <- "D:/L1polymORF/Data/BZ_NonRef"
NewL1RefOutPath       <- "D:/L1polymORF/Data/L1RefPacBioNA12878_DelRemoved.fa"
FastaFolder           <- "D:/L1polymORF/Data/"

# Set sequence piece length parameter (to identify reads reaching unique region)
SeqPcLength <- 30
MaxMisMatch <- 5

##############
# Load data
##############

# Load objects created by script AnalyzeReadsMapped2L1
load("D:/L1polymORF/Data/ReadsMapped2L1Info.RData")
L1ConsSeqDNASt <- DNAString(paste(L1ConsSeq, collapse = ""))
L1Length <- length(L1ConsSeqDNASt)
L1ConsLast500 <- subseq(L1ConsSeqDNASt, L1Length - 499, L1Length)

# Read reads mapped to the new L1 reference
ReadList <- scanBam("D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_sub_reads_aln2NewL1.sorted.minq1.bam")

# Create a boolean vector indicating which read is mapped to a new L1 insertion
blnNewL1 <- sapply(as.character(ReadList[[1]]$rname), function(x){
  SpltName <- strsplit(x, "_")[[1]]
  SpltName[length(SpltName)] == "New"
})

# Determine which reads should reach into the unique region right of 3' end
NrClipped <- sapply(ReadList[[1]]$cigar, NrClippedFromCigar)
bln3PEnd  <- NrClipped[2,] >= L1Length - FivePEnd
RL3PEnd   <- lapply(ReadList[[1]], function(y) y[bln3PEnd & blnNewL1])
NrClippedRight <- NrClipped[2, bln3PEnd & blnNewL1]
RL3PEnd$strand
RL3PEnd$seq[1:2]
c(RL3PEnd$seq[[1]], RL3PEnd$seq[[2]])

# Temporary exploration on strandedness of reads (TO BE REMOVED)
seq3P <- subseq(RL3PEnd$seq, width(RL3PEnd$seq) - NrClippedRight - 20, 
                width(RL3PEnd$seq) - NrClippedRight - 1)
PML <- lapply(1:length(seq3P), function(i) matchPattern(seq3P[[i]], 
          L1ConsSeqDNASt, max.mismatch = MaxMisMatch, with.indels=T))
bln1M <- 1* (sapply(PML, length) > 0)
seq3P[sapply(PML, length) > 0]
aggregate(bln1M ~ RL3PEnd$strand, FUN = mean)

# Loop through insertions that have reads reaching into unique region and test
# for each insertion whether part of the read maps onto unique region.
InsertionNames <- as.character(unique(RL3PEnd$rname))
InsertionWithMatches <- c()
PWAList <- c()
FastaFileNames <- c()
#InsertionName <- InsertionNames[1]
for (InsertionName in InsertionNames){
  
  cat("Analyzing insertion", InsertionName, "\n")
  SpltName <- strsplit(InsertionName, "_")[[1]]
  blnFullL1Info <- FullL1Info$chromosome == SpltName[1] &
    FullL1Info$start == as.numeric(SpltName[2])
  InsPos       <- FullL1Info$L1InsertionPosition.median[blnFullL1Info]
  TransdLength <- FullL1Info$L15PTransdSeq.median[blnFullL1Info]
  Strand       <- FullL1Info$L1Strand[blnFullL1Info]
  
  
  if (length(InsPos) > 0){
    if(!is.na(InsPos)){
      StartEnd  <- InsPos + c(0, (-2*(Strand == "+") + 1) * SeqPcLength)
      StartEnd500 <- InsPos + c(0, (-2*(Strand == "+") + 1) * 499)
      GR  <- GRanges(seqnames = SpltName[1], IRanges(start = min(StartEnd), 
                                                    end = max(StartEnd)))
      GR500 <- GRanges(seqnames = SpltName[1], IRanges(start = min(StartEnd500), 
                                                    end = max(StartEnd500)))
      SeqPiece  <- getSeq(BSgenome.Hsapiens.UCSC.hg19, GR)
      SeqPiece500 <- getSeq(BSgenome.Hsapiens.UCSC.hg19, GR500)
      SeqJunction <- c(L1ConsLast500, SeqPiece500[[1]])
      if (Strand == "-") {
        SeqPiece    <- reverseComplement(SeqPiece)
        SeqJunction <- c(L1ConsLast500, reverseComplement(SeqPiece500[[1]]))
      }
      blnInsertion <- RL3PEnd$rname == InsertionName
      cat("Analyzing", sum(blnInsertion), "reads\n")
      Seqs         <- RL3PEnd$seq[blnInsertion]
      SeqsNames    <- RL3PEnd$qname[blnInsertion]
      
      # Determine the position on the read that corresponds to the 3' end of 
      # the L1 (RightEnds)
      SeqLengths   <- width(Seqs)
      RightClips   <- NrClippedRight[blnInsertion] 
      RightEnds    <- SeqLengths - RightClips + L1Length - FivePEnd
      RightEnds500 <- RightEnds - 500
      blnREPos     <- RightEnds < SeqLengths
      cat(sum(blnREPos), "reads reaching into unique region\n")
      SeqEnds <- subseq(Seqs[blnREPos], RightEnds[blnREPos], SeqLengths[blnREPos]) 
      SeqEnds500 <- subseq(Seqs[blnREPos], RightEnds500[blnREPos], SeqLengths[blnREPos]) 
      PMList  <- lapply(SeqEnds, function(x) {
        matchPattern(SeqPiece[[1]], x, max.mismatch = MaxMisMatch, with.indels=T)
      })
      ML <- sapply(PMList, length)
      cat(sum(ML > 0), "reads match unique region\n")
      if (sum(ML) > 0){
        InsertionWithMatches <- c(InsertionWithMatches, InsertionName)
        LocalPWAList <- lapply(which(ML > 0), function(z){
          StartEnd <- InsPos + c(0, (-2*(Strand == "+") + 1) * width(SeqEnds)[z])
          GR <- GRanges(seqnames = SpltName[1], IRanges(start = min(StartEnd), 
                                                        end = max(StartEnd)))
          SeqPiece <- getSeq(BSgenome.Hsapiens.UCSC.hg19, GR)
          pairwiseAlignment(SeqEnds[z], SeqPiece, gapOpening = 4, gapExtension=2)
        })
        PWAList <- c(PWAList, LocalPWAList)
      }
      
      # Create a fasta file with the sequence of the junction and all reads 
      # spanning the junction
      SeqList <- list(RefSeq = s2c(as.character(SeqJunction)))
      for (j in 1:length(SeqEnds500)){
        ReadSeq <- SeqEnds500[[j]]
        SeqList <- c(SeqList, list(s2c(as.character(ReadSeq))))
      }
      names(SeqList) <- c("RefSeq", SeqsNames[blnREPos])
      FastaFileName <- paste(FastaFolder, InsertionName, ".fas", 
                                      sep = "")
      cat("Writing sequences to", FastaFileName, "\n")
      write.fasta(SeqList, names(SeqList), FastaFileName)
      FastaFileNames <- c(FastaFileNames, FastaFileName)
      
    }
  }
}
    
# Calculate the proportion of mismatching positions
sapply(1:length(PWAList), function(i) {
  length(PWAList[[i]]@pattern@mismatch[[1]]) / width(PWAList[[i]]@pattern@range)
})
sapply(1:length(PWAList), function(i) {
  length(PWAList[[i]]@pattern@indel[[1]]) / width(PWAList[[i]]@pattern@range)
})
sapply(1:length(PWAList), function(i) {
  width(PWAList[[i]]@subject@range)
})


# Loop through fasta files and align them
for (FastaFile in FastaFileNames){
  AlignedFile <- gsub(".fas", "_aligned.fas", FastaFile)
  cat("Aligning", FastaFile, "\n")
  run_MUSCLE(InputPath = FastaFile, OutputPath = AlignedFile) 
}

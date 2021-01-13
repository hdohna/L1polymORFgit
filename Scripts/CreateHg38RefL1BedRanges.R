# The following script creates a bed file for the non-polymorphic L1 ranges for the reference
# genome hg19 

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(GenomicRanges)
library(seqinr)

# Read in table with L1HS on Hg38
L1Table <- read.csv("D:/OneDrive - American University of Beirut/L1polymORF/Data/repeatsHg38_L1HS.csv",
                    as.is = T)
L1Table$ChrNr <- substr(L1Table$genoName, 4, nchar(L1Table$genoName))
L1Table$genoStartMinus1000 <- L1Table$genoStart - 1000
L1Table$genoEndPlus1000    <- L1Table$genoEnd + 1000
L1Table$idx <- 1:nrow(L1Table)
L1GR <- makeGRangesFromDataFrame(L1Table, seqnames.field = "genoName",
                                 start.field = "genoStart",
                                 end.field = "genoEnd")
blnFull <- width(L1GR) >= 6000
blnPlus <- as.vector(strand(L1GR)) == "+"
L1GR_left <- makeGRangesFromDataFrame(L1Table, 
                                      seqnames.field = "genoName",
                                      start.field = "genoStartMinus1000",
                                      end.field = "genoStart")

L1GR_right <- makeGRangesFromDataFrame(L1Table, 
                                       seqnames.field = "genoName",
                                       start.field = "genoEnd",
                                       end.field = "genoEndPlus1000")

export.bed(L1GR, "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HSRefRanges_hg38.bed")
export.bed(L1GR_left, "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HSRefRanges_hg38_left1000.bed")
export.bed(L1GR_right, "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HSRefRanges_hg38_right1000.bed")

# Read in L1 consensus 
L1Consens <- read.fasta("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1Hs_consensus.fas")
L1ConsensDNAST <- DNAString(paste(L1Consens[[1]], collapse = ""))

# Get L1HS sequences
L1HSSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, L1GR)

# Align L1HSSeq to consensus and calculate the proportion of nucleotides that differ from 
cat("Calculating the proportion of each L1 that differs from consensus ")
cat("(takes a few minutes) ...")
PropMismatch <- sapply(L1HSSeq, function(x){
  PwA <- pairwiseAlignment(x, L1ConsensDNAST)
  length(PwA@pattern@mismatch[[1]]) / length(x)
})
cat("done!\n")

# Get sequence of full sequences
L1FullGR  <- L1GR[blnFull]
L1FullSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, L1FullGR)
L1FullSeqList <- lapply(L1FullSeq, function(x) s2c(as.character(x)))
L1FullSeqList <- c(L1FullSeqList, list(L1Consens = toupper(L1FullSeqList[[1]])))

SeqNames <- paste(as.vector(seqnames(L1FullGR)), start(L1FullGR), 
                  end(L1FullGR), sep = "_")

# Write out full-length sequences and consensus sequence as fasta file
write.fasta(L1FullSeqList, names = c(SeqNames, "L1consensus"),
            file.out = "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1Hs_HG38_with_consensus_unaligned.fas")

# Save the objects
save.image(file = "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HSRefRanges_hg38.RData")
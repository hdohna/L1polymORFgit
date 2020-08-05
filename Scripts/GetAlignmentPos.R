# General description:
#    The following function maps genomic positions to positions within an alignment

# Arguments:
#   GrangesPos: genomic positions (as GRanges object)
#   GrangesL1: genomic ranges of L1 
#   AlignmentL1: alignment of L1s
# Output:
#   

#########################################################################
#                        TO BE DELETED:                                 #

source("file:///D:/L1polymORFgit/Scripts/_Start_L1polymORF.R")
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

library(rtracklayer)
library(ape)
library(seqinr)

# Read L1 catalog
L1Catalog <- read.csv("file:///D:/OneDrive - American University of Beirut/L1polymORF/Data/L1CatalogExtended.csv",
                      as.is = T)
L1CatalogSubset <- L1Catalog[!is.na(L1Catalog$start_HG38) & L1Catalog$Allele == 1 &
                             L1Catalog$blnInRef, ]

LiftoverList <- LiftoverL1Catalog(L1CatalogSubset,
                  ChainFilePath = "D:/OneDrive - American University of Beirut/L1polymORF/Data/hg38ToHg19.over.chain")
GRangesL1 <- LiftoverList$GRCatalogue_hg19
L1Seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, GRangesL1)

GR1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 10001, end = 10010))
GR2 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 10001, end = 10010), strand = "-")
getSeq(BSgenome.Hsapiens.UCSC.hg19, GR1)
getSeq(BSgenome.Hsapiens.UCSC.hg19, GR2)

# Read in L1 alignment
AlignmentL1 <- read.alignment("D:/OneDrive - American University of Beirut/L1Evolution/SawsanScriptsData/FinalScripts&RData/RData_SW/UpdatedFullSeq_withCon.fas",
                              format = "fasta")
AlignmentL1 <- read.fasta("D:/OneDrive - American University of Beirut/L1Evolution/SawsanScriptsData/FinalScripts&RData/RData_SW/UpdatedFullSeq_withCon.fas")
AlignNames <- attr(AlignmentL1,"name")
AlignAccAlleles <- sapply(AlignNames, function(x) {
  if (x == "L1Consensus") {
    c(Accession = "L1Consensus", Allele = "1")
  } else {
    Split <- strsplit(x, "\\.")[[1]]
    c(Accession = Split[1], Allele = Split[2])
    
  }
})
AlignmentL1NoAllele <- AlignmentL1[AlignAccAlleles["Allele", ] == "1"]
attr(AlignmentL1NoAllele,"name") <-  AlignAccAlleles["Accession", AlignAccAlleles["Allele", ] == "1"]

# Create a list of pairwise alignments
PairwiseAlignList <- lapply(1:length(L1Seqs), function(i){
  idxAlign <- which(attr(AlignmentL1NoAllele,"name") == GRangesL1$Acession[i])
  AlSeq <- AlignmentL1NoAllele[[idxAlign]]
  blnBlank <- AlSeq == "-"
  SeqStr <- paste(AlSeq[!blnBlank], collapse = "")
  pairwiseAlignment(L1Seqs[i], SeqStr)
})
as.vector(strand(GRangesL1)) == L1CatalogSubset$strand_L1toRef

# Repace alignments with negative scores by its reverse complement
PAScores <- sapply(PairwiseAlignList, function(PA) PA@score)
idxNeg   <- which(PAScores < 0)
for (i in idxNeg){
  idxAlign <- which(attr(AlignmentL1NoAllele,"name") == GRangesL1$Acession[i])
  AlSeq <- AlignmentL1NoAllele[[idxAlign]]
  blnBlank <- AlSeq == "-"
  SeqStr <- paste(AlSeq[!blnBlank], collapse = "")
  NewPA <- pairwiseAlignment(reverseComplement(L1Seqs[i]), SeqStr)
  PairwiseAlignList[[i]] <- NewPA
  cat("*************    ", i, "     *************\n")
  cat("Old strand:", L1CatalogSubset$strand_L1toRef[i], "\n")
  L1CatalogSubset$strand_L1toRef[i] <- switch(L1CatalogSubset$strand_L1toRef[i], `+` = '-', `-` = "+")
  cat("New strand:", L1CatalogSubset$strand_L1toRef[i], "\n")
}

# 
NIndelSubject <- sapply(PairwiseAlignList, function(PA) length(PA@subject@indel[[1]]))
NIndelPattern <- sapply(PairwiseAlignList, function(PA) length(PA@pattern@indel[[1]]))

which(NIndelSubject > 0 & NIndelPattern > 0)
hist(PAScores)
PA <- PairwiseAlignList[[4]]
PA@subject@indel
PA@subject@mismatch
PA@pattern@indel
PA@subject@unaligned[[1]][1:11]
PA@pattern@unaligned[[1]][1:11]
PA@subject@unaligned[[1]][6013:6023]
PA@pattern@unaligned[[1]][6013:6032]

OL <- findOverlaps(GRangesL1, GrangesPos)

##################################################################################

GetAlignmentPos <- function(Grange, Alignment, GenomeStarts, GenomeEnds) {
  
  if ((Grange@strand == "+")@values) { 
    SeqPos2 <- Grange@ranges@start + 1 - GenomeStarts
  } else {
    SeqPos2 <- 1 +GenomeEnds - Grange@ranges@start
  }
  
  # Row in data where chromosome, strand and range match with grange
  seqnb  <- which.min(abs(SeqPos2)) 
  PALign <- pairwiseAlignment(NewSeqHG19[[seqnb]], DNAString(as.character(data$L1Seq[[seqnb]])))
  patternStart <- PALign@pattern@range@start
  subjectStart <-  PALign@subject@range@start
  if(patternStart <subjectStart) # i.e indels in the beginning of pattern
  {SeqChr <- append(rep("-", subjectStart-patternStart),strsplit(as.character(PALign@pattern),"")[[1]])}
  else{SeqChr <- strsplit(as.character(PALign@pattern),"")[[1]]}
  SeqChr2 <- strsplit(as.character(PALign@subject),"")[[1]]
  
  SeqPos1 <- SeqPos2[seqnb]
  idxNuc2 <- which(SeqChr!="-")
  
  AlignPos2 <- idxNuc2[SeqPos1] #wrt to Palign
  
  SeqPos <- match(AlignPos2, which(SeqChr2 !="-")) # wrt to bac clone sequence in PAlign
  
  idxNuc <-  which(Alignment[seqnb,]!="-")
  
  AlignPos <- idxNuc[SeqPos]
  print(AlignPos)
  
}

# The following file reads an alignment of L1 sequences on hg19 and maps hg19 
# SNPs onto that alignment

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Source required packages
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(seqinr)


# Sequence of start and end of ORF1 and ORF2
StartSeqORF1 <- "ATGGGGAAA"
StartSeqORF2 <- "ATGACAGGA"
EndSeqORF1   <- "GCCAAAATGTAA"
EndSeqORF2   <- "GGTGGGAATTGA"

# Read repeat masker table for L1HS
L1Table <- read.csv("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HS_repeat_table_Hg19.csv", as.is = T)

# Create GRanges objects with L1 Seqences
L1GR <- makeGRangesFromDataFrame(L1Table, seqnames.field = "genoName",
                                 start.field = "genoStart",
                                 end.field = "genoEnd")

# Get sequences and create a character of sequence names
L1Seq    <- getSeq(BSgenome.Hsapiens.UCSC.hg19, L1GR)
SeqNames <- paste(as.vector(seqnames(L1GR)), start(L1GR), end(L1GR), sep = "_")

# Form different subsets and write them out as fasta files
L1Aligned <- read.fasta(file = "file:///D:/OneDrive - American University of Beirut/L1polymORF/Data/L1seqHg19_minLength6000_aligned.txt")
# L1Aligned <- read.alignment(file = "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1seqHg19_minLength6000_aligned.txt",
#                             format = "fasta")

# Read vcf with variants in LINE-1s 
L1Variants  <- ReadVCF("D:/OneDrive - American University of Beirut/L1polymORF/Data/VariantsInL1.recode.vcf")
L1Variants$chromosome <- paste("chr", L1Variants$X.CHROM, sep = "")

# Create a GRanges object of variants inside L1s and their flanking regions
L1VarGR <- makeGRangesFromDataFrame(L1Variants, 
                                    start.field = "POS",
                                    end.field = "POS")

# Get indicies of sequenced positions
idxSeqPosList <- lapply(L1Aligned, function(x) which(x != "-"))

# Identify start and end of ORF1 and ORF2 on alignment
i <- 1
ORFStartEndList <- lapply(1:length(L1Aligned), function(i){
  Seq <- L1Aligned[[i]]
  idxSeq <-  which(Seq != "-")
  ORF1StartMatch <- matchPattern(StartSeqORF1, paste(toupper(Seq[idxSeq]), collapse = ""))
  ORF1EndMatch   <- matchPattern(EndSeqORF1, paste(toupper(Seq[idxSeq]), collapse = ""))
  ORF2StartMatch <- matchPattern(StartSeqORF2, paste(toupper(Seq[idxSeq]), collapse = ""))
  ORF2EndMatch   <- matchPattern(EndSeqORF2, paste(toupper(Seq[idxSeq]), collapse = ""))
  list(ORF1Start = idxSeq[start(ORF1StartMatch)],
       ORF1End = idxSeq[end(ORF1EndMatch)],
       ORF2Start = idxSeq[start(ORF2StartMatch)],
       ORF2End = idxSeq[end(ORF2EndMatch)])
  
})
table(unlist(lapply(ORFStartEndList, function(x) x$ORF1Start)))
table(unlist(lapply(ORFStartEndList, function(x) x$ORF1End)))
table(unlist(lapply(ORFStartEndList, function(x) x$ORF2Start)))
table(unlist(lapply(ORFStartEndList, function(x) x$ORF2End)))

# Determine overlaps between SNPs and LINE-1
OL <- findOverlaps(L1GR, L1VarGR)
idxMinus <- 1 + c(as.vector(strand(L1GR)) == "-")
StartEnd <- cbind(start(L1GR), end(L1GR))
L1Starts <- sapply(seq_along(idxMinus), function(i) StartEnd[i, idxMinus[i]])

SNPpos <- abs(start(L1Starts)[OL@from] - start(L1VarGR)[OL@to]) + 1


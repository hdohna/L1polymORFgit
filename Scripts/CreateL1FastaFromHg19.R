# The following file reads in a repeatmasker table of L1s and creates a fasta file
# with full-length L1 

# Source required packages
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(seqinr)

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
bln6000 <- width(L1Seq) >= 6000
L1List <-  lapply(which(bln6000), function(x) s2c(as.character(L1Seq[x])))
write.fasta(L1List, names = SeqNames[bln6000], 
            file.out = "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1seqHg19_minLength6000_unaligned.fas")
bln500 <- width(L1Seq) >= 500
L1List <-  lapply(which(bln500), function(x) s2c(as.character(L1Seq[x])))
write.fasta(L1List, names = SeqNames[bln500], 
            file.out = "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1seqHg19_minLength500_unaligned.fas")


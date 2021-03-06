##############################################
#
# General description:
#
#   The following script reads in data sources genomic regions , r

# Input:
#
#    D:/L1polymORF/Data/repeatsHg38: table with all repeats
#   

# Output:
#   
#    L1HS_repeat_table.csv: csv file with all L1HS repeats
#    L1Sequences_reptab.fas: fasta file with all L1 sequences

##############################################

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.r')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg19)
library(seqinr)

# Load consensus sequence
L1Consens <- read.fasta("D:/OneDrive - American University of Beirut/L1polymORF/Data/Homo_sapiens_L1_consensus.fa")
L1ConsensDNAST <- DNAString(paste(L1Consens[[1]], collapse = ""))

# Read in repeatMasker table 
L1Table   <- read.csv("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HS_repeat_table_Hg19.csv")

# Genomic ranges for each dataset
L1GR <- makeGRangesFromDataFrame(L1Table, 
                                                   seqnames.field = "genoName", 
                                                   start.field ="genoStart", 
                                                   end.field = "genoEnd", 
                                                   strand.field = "strand")

# Get L1HS sequences
L1HSSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, L1GR)

# Align L1HSSeq to consensus and calculate the proportion of nucleotides that differ from 
cat("Calculating the proportion of each L1 that differs from consensus ")
cat("(takes a few minutes) ...")
PropMismatch <- sapply(L1HSSeq, function(x){
  PwA <- pairwiseAlignment(x, L1ConsensDNAST)
  length(PwA@pattern@mismatch[[1]]) / length(x)
})
cat("done!\n")

save.image(file = "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HS_PropMismatch.RData")

# Test for correlation between proportion of mismatch and length of the L1
cor.test(width(L1HSSeq), PropMismatch, method = "spearman")
blnFull <- width(L1HSSeq) > 6000
blnAbove2000 <- width(L1HSSeq) > 2000
cor.test(width(L1HSSeq)[!blnFull], PropMismatch[!blnFull], method = "pearson")
t.test(PropMismatch ~ blnFull)
wilcox.test(PropMismatch ~ blnFull)
wilcox.test(PropMismatch[!blnFull] ~ blnAbove2000[!blnFull])
wilcox.test(PropMismatch[blnAbove2000] ~ blnFull[blnAbove2000])
boxplot(PropMismatch ~ blnFull,ylim = c(0, 0.03))
sum(!blnFull)
sum((!blnFull) & blnAbove2000)

# Turn sequences them into a list and write it out as fasta file
L1HSSeqList <- lapply(L1HSSeq, function(x) tolower(s2c(as.character(x))))
names(L1HSSeqList) <- paste(L1Table$genoName, L1Table$genoStart,
                            L1Table$genoEnd, L1Table$strand, sep = "_")
write.fasta(c(L1Consens, L1HSSeqList), c("L1consensus", names(L1HSSeqList)),
            file.out = "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HS_withConsens_unaligned.fas")
            
write.fasta(list(L1Consens[[1]], tolower(s2c(as.character(L1HSSeq[1])))), 
                 c("L1consensus", "L1fragm"),
            file.out = "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HS_example.fas")

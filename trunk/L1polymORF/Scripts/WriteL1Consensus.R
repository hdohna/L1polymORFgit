# The following script reads L1 sequences (obtained from http://www.girinst.org/repbase/)
# and selects one as the consensus sequence
library(ape)
library(seqinr)
L1seqs <- read.dna("D:/L1polymORF/Data/Homo_sapiens_L1", format = "fasta",
                   as.character = T)
#L1seqs <- read.fasta("D:/L1polymORF/Data/Homo_sapiens_L1")
names(L1seqs) <- gsub("\t", "_", names(L1seqs))
names(L1seqs) <- gsub(" ", "_", names(L1seqs))
sapply(L1seqs, length)
write.fasta(L1seqs$L1HS_L1_Homo_sapiens, "L1HS_L1_Homo_sapiens",
            "D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fa")

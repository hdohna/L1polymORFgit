# The following script reads L1 sequences (obtained from http://www.girinst.org/repbase/)
# and selects the L1Hs consensus sequence
library(ape)
library(seqinr)
L1seqs <- read.dna("D:/L1polymORF/Data/Homo_sapiens_L1", format = "fasta",
                   as.character = T)

# Remove \t from names
names(L1seqs) <- gsub("\t", "_", names(L1seqs))
names(L1seqs) <- gsub(" ", "_", names(L1seqs))
SeqLengths <- sapply(L1seqs, length)
hist(SeqLengths)
any(duplicated(names(L1seqs)))

# Write out L1HS
write.fasta(L1seqs$L1HS_L1_Homo_sapiens, "L1HS_L1_Homo_sapiens",
            "D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fa")

# Find L1PA with maximum length
idxL1HS <- which(names(L1seqs) == "L1HS_L1_Homo_sapiens")
idxL1PA <- grep("L1PA", names(L1seqs))
idxL1PA_O  <- order(SeqLengths[idxL1PA], decreasing = T)
idxMaxL1PA  <- idxL1PA[which.max(SeqLengths[idxL1PA])]
idxHighL1PA <- idxL1PA[idxL1PA_O[1:5]]

# Write out L1HS and L1PA
write.fasta(L1seqs[c(idxL1HS, idxMaxL1PA)], 
            names(L1seqs)[c(idxL1HS, idxMaxL1PA)],
            "D:/L1polymORF/Data/Homo_sapiens_L1HSL1PA_consensus.fa")
write.fasta(L1seqs[c(idxL1HS, idxHighL1PA)], 
            names(L1seqs)[c(idxL1HS, idxHighL1PA)],
            "D:/L1polymORF/Data/Homo_sapiens_L1HSL1PA_Multi.fa")

write.fasta(L1seqs[c(idxL1HS, idxHighL1PA)], 
            names(L1seqs)[c(idxL1HS, idxHighL1PA)],
            "D:/L1polymORF/Data/Homo_sapiens_L1HSL1PA_Multi.fa")

##############################################
#
# General description:
#
#   The following script reads in accession numbers of clones that contain 
#   full-length L1 insertions and gets L1 sequences from the clones

# Input:
#
#    D:/OneDrive - American University of Beirut/L1polymORF/Data/L1CatalogExtended.csv: 
#    table with all full-length L1
#   

# Output:
#   
#    L1HSFromClones,fas

##############################################

########################################
#                                      #
#  Source packages and set parameters  #
#                                      #
########################################

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.r')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(seqinr)

#######################################
#                                     #
#     Read data                       #
#                                     #
#######################################

# Load consensus sequence
L1Consens <- read.fasta("D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fa")
L1Consens <- paste(L1Consens[[1]], collapse = "")
L1Consens <- DNAString(L1Consens)

# Read in catalog
L1Catalog <- read.csv("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1CatalogExtended.csv",
                      as.is = T)
DuplAcc <- duplicated(L1Catalog$Accession)

# Read in clone sequences
CloneSeqs <- read.fasta("D:/L1polymORF/Data/L1CloneSeq.fas")
CloneSeqChar <- sapply(CloneSeqs, function(x) paste(x, collapse = ""))
names(CloneSeqs)
AccMatch <-  match(names(CloneSeqs), L1Catalog$Accession)
NewDF <- getNonRefL1(L1Consens = L1Consens,
                    AccNrs = names(CloneSeqs), 
                    Seqs = DNAStringSet(CloneSeqChar), 
                    MinMatchWidth = 5500, 
                    FlankSize = 100,
                    Chromosomes = L1Catalog$Chromosome[AccMatch], 
                    RefGenome = BSgenome.Hsapiens.UCSC.hg38,
                    blnLocateL1inRef = F)

# Write the L1 sequences out as a fasta file 
blnNotNA <- ! is.na(NewDF$L1Seq)
Seqs <- lapply(c(as.character(L1Consens), NewDF$L1Seq[blnNotNA]), function(x){
  strsplit(x, "")[[1]]
})
write.fasta(sequences = Seqs, names = c("L1Consensus", names(CloneSeqs)[blnNotNA]), 
            file.out = "D:/OneDrive - American University of Beirut/L1Evolution/Data/L1SeqsFromClonesUnaligned.fas")

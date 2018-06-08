# The following script reads in a version of the L1 catalog, adds columns
# and saves the result as "L1CatalogExtended.csv:

# Load packages
library(ape)
library(seqinr)

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Read in table with known L1 
L1Catalog <- read.csv("D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv", as.is = T)

##########################
#                        #
#      Add columns       #
#                        #
##########################

# Add numeric activity
L1Catalog$ActivityNum <- L1Catalog$Activity
L1Catalog$ActivityNum <- gsub("<", "", L1Catalog$ActivityNum)
L1Catalog$ActivityNum <- as.numeric(L1Catalog$ActivityNum)

# Add numeric allele frequency
L1Catalog$Allele_frequency_Num <- as.numeric(L1Catalog$Allele_frequency)

# Add boolean variable indicating L1s in reference genopme
L1Catalog$blnInRef <- (L1Catalog$end_HG38 - L1Catalog$start_HG38) > 6000 

# Write out extended catalog
write.csv(L1Catalog, "D:/L1polymORF/Data/L1CatalogExtended.csv")


#################################
#                               #
#     Write out sequences       #
#                               #
#################################

# Read prototypic L1 sequences (from http://www.girinst.org/repbase/)  
L1consens <- read.dna("D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fa", format = "fasta",
                   as.character = T)

# Put L1 sequences from catalog in list
idxSeq <- which(!is.na(L1Catalog$L1Seq))
L1HS_SeqList <- lapply(L1Catalog$L1Seq[idxSeq],
                       function(x) tolower(s2c(x)))

# Prepend consensus sequence
L1HS_withConsens <- c(list(L1consens[1,]), L1HS_SeqList)
L1CatNames <- paste(L1Catalog$Accession[idxSeq], 
                    L1Catalog$Allele[idxSeq], sep = "_")
names(L1HS_withConsens) <- c("L1consensus", L1CatNames)

which(sapply(L1HS_withConsens, function(x) any(x == "I")))

# Write sequences as fasta file and align
write.fasta(L1HS_withConsens, names(L1HS_withConsens), 
            file.out = "D:/L1polymORF/Data/L1Catalog_Unaligned_withConsens.fas")
# run_MUSCLE(
#   InputPath = "D:/L1polymORF/Data/L1Catalog_Unaligned_withConsens.fas", 
#            OutputPath = "D:/L1polymORF/Data/L1Catalog_Aligned_withConsens.fas")

# Read alignment and write it out as nexus file
L1HSAligned <- read.dna("D:/L1polymORF/Data/L1Catalog_Aligned_withConsens.fas",
                        format = "fasta", as.character = T)

# Get indices of consensus sequence
idxConsens <- which(L1HSAligned$L1consensus != "-")
length(idxConsens)
idxCatSeq  <- which(names(L1HSAligned) != "L1consensus")
NamesCatSeq  <- names(L1HSAligned)[idxCatSeq]
CatMatch   <- match(NamesCatSeq, paste(L1Catalog$Accession, 
                                                L1Catalog$Allele, sep = "_"))

# Get start and end of each sequence with respect to the reference sequence
idxStartEndSeq <- sapply(idxCatSeq, function(x) {
  idxSeq <- intersect(idxConsens, which(L1HSAligned[[x]] != "-"))
  c(Start = which(idxConsens == min(idxSeq)), 
    End = which(idxConsens == max(idxSeq)))
})

# Get the highest start and lowest end of all sequences with non-zero transposition
# activity
ActNum <- L1Catalog$ActivityNum[CatMatch]
idxActNonzero <- which(ActNum > 1)
MaxStart <- max(idxStartEndSeq["Start", idxActNonzero])
MinEnd   <- min(idxStartEndSeq["End", idxActNonzero])

plot(jitter(idxStartEndSeq["Start", ]), ActNum)
plot(jitter(idxStartEndSeq["End", ]), ActNum)

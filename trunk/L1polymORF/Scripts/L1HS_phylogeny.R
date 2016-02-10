##############################################
#
# General description:
#
#   The following script reads an alignment of L1HS sequnces and constructs a
#   phylogeny 

# Input:
#
#    L1Sequences_reptab_alignedMAFF.fas: fasta file with all L1 sequences
#   

# Output:
#   
#    : 

##############################################

########################################
#                                      #
#  Source packages and set parameters  #
#                                      #
########################################

# Source start script
#source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(seqinr)
library(ape)

# Files and folders
AlignFileName <- "D:/L1polymORF/Data/L1Sequences_reptab_alignedMAFF"

# Source all functions from Functions folder
AllFunctions <- list.files(path = "D:/GeneralRFunctions/", 
                           pattern = ".[rR]", full.names = T)
sapply(AllFunctions, source)

#######################################
#                                     #
#     Read bed file and               #
#   save fasta file with sequences    #
#                                     #
#######################################

# Read in alignment
L1HSAlign <- read.dna(AlignFileName, format = "fasta")

# Create a tree and plot it
Tree <- bionj(dist.dna(L1HSAlign))
plot(Tree, cex = 0.075, edge.width = 0.1)
CreateDisplayPdf("D:/L1polymORF/Figures/L1Tree.pdf")

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
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(seqinr)
library(ape)

#######################################
#                                     #
#   Read tree and activity values     #
#                                     #
#######################################

# Read activity table
ActTable <- read.table("D:/L1polymORF/Data/ActivityBT_L1Catalogue.txt",
                       col.names = c("Name", "Activity"))

# Read in MrBayes consensus tree
Alignment <- read.dna("D:/L1polymORF/Data/L1Catalogue_Aligned_withRoot.fas",
                      format = "fasta")

# Root tree by L1PA sequence and remove it
Tree <- bionj(dist.dna(Alignment))
Tree <- root(Tree, which(Tree$tip.label == "L1PA"))
Tree <- drop.tip(Tree, which(Tree$tip.label == "L1PA"))

# Color tree by activity
ActMatch <- match(Tree$tip.label, ActTable$Name)
Act      <- ActTable$Activity[ActMatch]
Edgecolors <- rep("black", length(Tree$edge.length))

# Get all indices of tree edges that lead to tips with a particular subtype
LowEdges <- which.edge(Tree, which(Act  == 1)) 
MediumEdges <- which.edge(Tree, which(Act  == 2)) 
HighEdges   <- which.edge(Tree, which(Act  == 3)) 

# Color edges from H1 subtypes red and from H2 blue
Edgecolors[LowEdges] <- "yellow"
Edgecolors[MediumEdges] <- "orange"
Edgecolors[HighEdges]   <- "red"

plot(Tree, cex = 0.05, edge.width = 0.1, edge.color = Edgecolors)
CreateDisplayPdf("D:/L1polymORF/Figures/L1TreeAct_Catalogue.pdf")

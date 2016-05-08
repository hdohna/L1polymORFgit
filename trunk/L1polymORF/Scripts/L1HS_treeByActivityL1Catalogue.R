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

############################################
#                                          #
#   Read alignment and activity values     #
#                                          #
############################################

# Read activity table
ActTable <- read.table("D:/L1polymORF/Data/ActivityBT_L1Catalogue.txt",
                       col.names = c("Name", "Activity"))

# Read in MrBayes consensus tree
Alignment <- read.dna("D:/L1polymORF/Data/L1Catalogue_Aligned_withRoot.fas",
                      format = "fasta")

# Read in table with known L1 
L1Catalogue <- read.csv("D:/L1polymORF/Data/L1Catalogue_Fri_Apr_01_18-28-08_2016.csv", 
                        as.is = T)

############################################
#                                          #
#  Construct tree and plot activity values #
#                                          #
############################################

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

#################################################
#                                               #
#  Compare activity values to closest distance  #
#                                               #
#################################################

# Create a matrix of pairwise proportion of nucleotide differences
DNADist <- dist.dna(Alignment, as.matrix = T)
diag(DNADist) <- NA

# Get the allele of each matrix element
DistAllele <- sapply(row.names(DNADist), function(x) strsplit(x, "_")[[1]][2])

# Subset distance matrix to get only allele 1 per insertion and to remove outlier 
DNADist <- DNADist[DistAllele == "1", DistAllele == "1"]
DNADist <- DNADist[rownames(DNADist) != "AC117496_1", 
                   rownames(DNADist) != "AC117496_1"]

# Get accession numbers and match activities by accession numbers
DistAcc  <- sapply(row.names(DNADist), function(x) strsplit(x, "_")[[1]][1])
ActMatch <- match(DistAcc, L1Catalogue$Accession)
Act      <- L1Catalogue$Activity[ActMatch]

# Calculate for each L1 its distance to the closest other L1
Mindist <- sapply(1:nrow(DNADist), function(x){
  min(DNADist[x,], na.rm = T)
})

# Plot activity vs distance and test for correlations
Act[Act == "<1"] <- "0.5"
ActNum <- as.numeric(Act) 
plot(ActNum, Mindist, xlab = "Activity", ylab = "Nuc diff. to closest L1")
plot(ActNum, Mindist, xlab = "Activity", ylab = "Nuc diff. to closest L1",
     ylim = c(0, 0.02))
text(150, 0.015, "P < 0.0001")
CreateDisplayPdf("D:/L1polymORF/Figures/L1ActivityVsMinDist.pdf")
cor(ActNum, Mindist, use = "pairwise.complete.obs")
cor.test(ActNum, Mindist, use = "pairwise.complete.obs")

# Test for correlation when all zero activities are removed
cor.test(ActNum[ActNum > 0], Mindist[ActNum > 0], 
         use = "pairwise.complete.obs")



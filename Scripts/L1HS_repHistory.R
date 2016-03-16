##############################################
#
# General description:
#
#   The following script reads an alignment of L1HS sequences (full-length and 
#   fragments) and assigns each fragment to its closest full-length

# Input:
#
#    L1HS_L100_aligned.fas: fasta file with all L1 sequences
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

# Files and folders
AlignFileName <- "D:/L1polymORF/Data/L1HS_L100_aligned.fas"

# Minimum number of associated fragments to be included in analysis
MinFrags <- 0

#######################################
#                                     #
#     Read bed file and               #
#   save fasta file with sequences    #
#                                     #
#######################################

# Read in alignment
L1HSAlign <- read.dna(AlignFileName, format = "fasta")

# Create a tree and plot it
DistL1 <- dist.dna(L1HSAlign, model = "raw", pairwise.deletion = T)
sum(is.na(DistL1))
sum(!is.na(DistL1))

# Count nr nucs and get all indices of full-length 
L1HSAlignMat <- as.character(L1HSAlign)
NrNuc <- apply(L1HSAlignMat, 1, FUN = function(x) sum(x != "-"))
blnFullLength <- NrNuc > 6000

# # Create a plot of nucleotides ordered by framgment length
# NrNucPerCol <- apply(L1HSAlignMat, 2, FUN = function(x) sum(x != "-"))
# L1AlMatReduced <- L1HSAlignMat[,NrNucPerCol > 100]
# NrNuc <- apply(L1AlMatReduced, 1, FUN = function(x) sum(x != "-"))
# NrNucRank <- nrow(L1AlMatReduced) - rank(NrNuc, ties.method = "first")
# 
# 
# plot(c(1, ncol(L1AlMatReduced)), c(0, nrow(L1AlMatReduced)), 
#      type = "n", xaxt = "n", yaxt = "n", xlab = "",
#      ylab = "", bty = "n")
# 
# for (i in 1:nrow(L1AlMatReduced)){
#   XVals <- which(L1AlMatReduced[i,]!= "-")
# #   lines(c(1, ncol(L1AlMatReduced)), c(NrNucRank[i], NrNucRank[i]),
# #         col = "grey", lwd = 0.2)
#   points(XVals, rep(NrNucRank[i], length(XVals)), pch = 15,
#          cex = 0.1, col = grey(0.5, 0.5))
# }
# lines(colSums(L1AlMatReduced != "-"), col = "red",
#       lwd = 2)
# CreateDisplayPdf('D:/L1polymORF/Figures/AlignmentPlot.pdf')

# Get for each fragment its closest full-length
idxFullLength <- which(blnFullLength)
Fragmnet2FullLength <- t(sapply(which(!blnFullLength), function(i){
  Dist2Full <- DistL1[PtIndices2DistVect(i, idxFullLength, nrow(L1HSAlignMat))]
  idxFull <- idxFullLength[which.min(Dist2Full)]
  names(idxFull) <- NULL
  c(idxFrgm = i, idxFull = idxFull)
}))

# Get for each full-length the distance to all the closest framgments
FullDist2Frag <- lapply(idxFullLength, function(i){
  CurrentSubset  <- Fragmnet2FullLength[,"idxFull"] == i
  FragSubset     <- Fragmnet2FullLength[CurrentSubset,"idxFrgm"]
  DistL1[PtIndices2DistVect(i, FragSubset, nrow(L1HSAlignMat))]
})
names(FullDist2Frag) <- names(NrNuc)[idxFullLength]

NrFrags <- sapply(FullDist2Frag, length)
max(NrFrags)
max(unlist(FullDist2Frag))

# Create distance classes
DiffClasses <- seq(0, 0.52, 0.002)
idxSomeFrags <- which(NrFrags >= MinFrags)
Cols <- rainbow(length(idxSomeFrags))
plot(DiffClasses, DiffClasses, ylim = c(0, 15), type = "n",
     xlim = c(0, 0.05), xlab = "Difference to full-length",
     ylab = "Count")

# Create a matrix that counts for each full-length L1 that has some fragments 
# associated the number of associated fragments per distance class
CountMat <- sapply(1:length(idxSomeFrags), function(i){
  x <- FullDist2Frag[[idxSomeFrags[i]]]
  H <- hist(x, breaks = DiffClasses, plot = F)
  lines(H$mids, H$counts, col = Cols[i])
  H$counts
})
colnames(CountMat) <- names(NrFrags)[idxSomeFrags]
rownames(CountMat) <- DiffClasses[-1]
CreateDisplayPdf('D:/L1polymORF/Figures/L1repHistories.pdf')

# Save results
save.image(file = "D:/L1polymORF/Data/L1repHistoryResults.RData")


# The script below reads GRanges with L1s from the 1000 Genome data
# Phase 3, available at ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
# The script takes matches L1 insertions to 
# The Granges were created by the script 'Create_1000G_L1GRanges.R'

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(DescTools)
library(ade4)
library(vegan)
library(ggplot2)
library(grid)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Specify file paths
G1000SamplePath <- 'D:/L1polymORF/Data/1000GenomeSampleInfo.txt'
L1GRPath        <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
L1RefRangePath  <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'

# Specify parameters
NrInfoCols <- 9

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

# Load previously generated objects
load(L1RefRangePath)
load(L1GRPath)

# Read in table with Info about 1000 genome samples 
SampleInfo_1000Genome <- read.table(G1000SamplePath, as.is = T, header = T)

# # Get a list of samples per population
SamplePerPopList <-  lapply(unique(SampleInfo_1000Genome$super_pop), function(Pop){
  blnPop         <- SampleInfo_1000Genome$super_pop == Pop
  which(SampleColumns %in% SampleInfo_1000Genome$sample[blnPop])
})

PCA <- princomp(t(L1_1000G[1:1000,SampleColumns]))
plot(PCA$sdev[1:10]/sum(PCA$sdev), xlab = "Rank")
mp <- barplot(PCA$sdev[1:20]^2/sum(PCA$sdev^2), xlab = "Component rank",
        ylab = "Proportion of variance", xaxt = "n")
axis(1, at = mp[c(1, 5, 10, 15, 20)], c(1, 5, 10, 15, 20))
dim(PCA$loadings)
dim(PCA$scores)
plot(PCA$scores[,1], PCA$scores[,2])
Cols <- rainbow(length(SamplePerPopList))
for (i in 1:length(SamplePerPopList)){
  idx <- SamplePerPopList[[i]]
  points(PCA$scores[idx, 1], PCA$scores[idx, 2], col = Cols[i])
}

plot(PCA)

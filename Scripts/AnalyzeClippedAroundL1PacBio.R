# The following script test for each L1 in the catalog whether it is present
# in a genome using PacBio data

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(csaw)
library(ShortRead)
library(Rsamtools)

# Load data produced by GetNrClippedAroundL1PacBio
load(file = "D:/L1polymORF/Data/ClippedReadsL1_PacBio.RData")

# Get mean number clipped left and right
Left5P <- sapply(ClipList, function(x) x$MeanClippedLeft["CigNr_5Prime",])
Left3P <- sapply(ClipList, function(x) x$MeanClippedLeft["CigNr_3Prime",])
Right5P <- sapply(ClipList, function(x) x$MeanClippedRight["CigNr_5Prime",])
Right3P <- sapply(ClipList, function(x) x$MeanClippedRight["CigNr_3Prime",])

par(mfrow = c(2, 2))
plot(Left5P[1,], type = "l", ylim = c(min(Left5P), max(Left5P)))
for (i in 1:nrow(Left5P)){
  lines(Left5P[i,])
}
plot(Left3P[1,], type = "l", ylim = c(min(Left3P), max(Left3P)))
for (i in 1:nrow(Left5P)){
  lines(Left3P[i,])
}
plot(Right5P[1,], type = "l", ylim = c(min(Right5P), max(Right5P)))
for (i in 1:nrow(Left5P)){
  lines(Right5P[i,])
}
plot(Right3P[1,], type = "l", ylim = c(min(Right3P), max(Right3P)))
for (i in 1:nrow(Right3P)){
  lines(Right3P[i,])
}

RatioL5P <- Left5P[,1]/rowMeans(Left5P[,-1])
RatioL3P <- Left3P[,1]/rowMeans(Left3P[,-1])
RatioR5P <- Right5P[,1]/rowMeans(Right5P[,-1])
RatioR3P <- Right3P[,1]/rowMeans(Right3P[,-1])

par(mfrow = c(2, 2))
hist(RatioL5P)
hist(RatioL3P)
hist(RatioR5P)
hist(RatioR3P)

par(mfrow = c(2, 2))
plot(RatioL5P, RatioR3P)
plot(RatioL5P, RatioL3P)
plot(RatioL3P, RatioR5P)
plot(RatioR5P, RatioR3P)

cor(RatioL5P, RatioR3P)
cor(RatioL5P, RatioL3P)
cor(RatioL3P, RatioR5P)
cor(RatioR5P, RatioR3P)


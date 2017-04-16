# The following script determines how the number of detected L1 insertions 
# changes with the number of samples

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')
NrInfoCols <- 9

# Load 100 genome data
L1Table <- read.table("D:/L1polymORF/Data/L1_1000G_withGenoNum")

# Subset to get full-length insertions
L1InfoFull    <- L1Table[L1Table$InsLength > 5950, ]
SampleColumns <- colnames(L1Table)[(NrInfoCols + 1):ncol(L1Table)]

# Create a list that contains for each individual the L1s he/she has
L1List <- lapply(SampleColumns, function(x) which(L1InfoFull[,x] > 0))

# Create a matrix where the column corresponds to the number of individuals 
# sampled, the row is a random iteration and the entry is the number of
# distinct L1s
NrIter      <- 1000
NrPerSample <- 50
SampleMat   <- matrix(nrow = NrIter, ncol = length(SampleColumns) %/% NrPerSample)
ListIndices <- 1:length(L1List)
SampleSeq   <- 1:ncol(SampleMat) * NrPerSample
for (i in 1:NrIter){
  RandomSample <- sample(ListIndices, length(ListIndices), replace = T)
  SampleMat[i,] <- sapply(SampleSeq, function(x) {
    length(unique(unlist(L1List[RandomSample[1:x]])))})
}

# Plot median and quantiles
Cols <- rainbow(10)
QuantileMat <- apply(SampleMat, 2, FUN = function(x) quantile(x, c(0.05, 0.5, 0.95)))
idxFw <- 1:ncol(SampleMat)
idxRv <- ncol(SampleMat):1
plot(QuantileMat[2,], type = "n", ylim = c(150, max(QuantileMat)), 
     ylab = 'Number of L1 insertions detected', xaxt = "n", 
     xlab = "Sample size")
TickPos <- seq(1, length(SampleSeq), 10) - 1
TickPos <- TickPos[-1]
axis(1, at = TickPos, labels = SampleSeq[TickPos])
polygon(c(idxFw, idxRv), c(QuantileMat[1, idxFw], QuantileMat[3, idxRv]),
        col = Cols[7], border = NA)
#lines(QuantileMat[2,], lwd = 2, col = "red")
lines(colMeans(SampleMat), lwd = 2, col = "red")
CreateDisplayPdf("D:/L1polymORF/Figures/L1VsNumberSampled.pdf",
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7, width = 7)

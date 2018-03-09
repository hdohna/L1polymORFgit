# Load in file with 1000 genome L1 singletons
SingleL1 <- read.delim("D:/L1polymORF/Data/Singleton_L1_chr1", sep = " ",
                       header = F)

# 1000 genome data
load("D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData")

# Getting all singletons from L1_1000G
blnSingle <- rowSums(L1_1000G[,SampleColumns]) == 1 & L1_1000G$CHROM == 1
sum(blnSingle)

idxCol <- sapply(which(blnSingle), function(x){
  which(L1_1000G[x,SampleColumns] == 1)
})

PredictCol <- ((SingleL1$V10 - sum(nchar(SampleColumns)) - 2503) %% (4*2504))/4
PredictCol <- ceiling((SingleL1$V10 %% 10016)/4)

idxCol[idxCol != PredictCol]
idxCol - PredictCol
plot(idxCol, PredictCol)
lines(c(0, 10^10), c(0, 10^10))

L1_1000G_reduced[which(blnSingle)[idxCol != ceiling(PredictCol)],]
L1_1000G_reduced[blnSingle,]
 plot(which(blnSingle))
2158 / 346

for (i in 4){
  plot((SingleL1$V10 %% (i*2504))/i, idxCol, main = i)
  lines(c(0, 10^10), c(0, 10^10))
}


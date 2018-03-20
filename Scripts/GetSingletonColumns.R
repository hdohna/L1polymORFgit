# 1000 genome data
load("D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData")

# Specify chromosome
Chr <- 12

# Load in file with 1000 genome L1 singletons
SingleL1 <- read.delim(paste("D:/L1polymORF/Data/Singleton_L1_chr", Chr, sep = ""), 
                       sep = " ", header = F)


# Getting all singletons from L1_1000G
blnSingle <- rowSums(L1_1000G[,SampleColumns]) == 1 & L1_1000G$CHROM == Chr
sum(blnSingle)

# Get index of singleton column 
idxCol <- sapply(which(blnSingle), function(x){
  which(L1_1000G[x,SampleColumns] == 1)
})

# # Predict column from singleton file
# minusFirstRow <- SingleL1$V10 - 20031
# beyondRow     <- minusFirstRow %/% 10016 >= 5226911
# PredictCol    <- ceiling(((minusFirstRow - beyondRow*11401) %% 10016)/4)
# 
# idxCol[idxCol != PredictCol]
# idxCol - PredictCol
# plot(idxCol, PredictCol)
# lines(c(0, 10^10), c(0, 10^10))

ColPredict <- GetSingletonColumns(SingleL1)
plot(idxCol, ColPredict$Col)
lines(c(0, 10^10), c(0, 10^10))
table(ColPredict$Allele)

max(abs(idxCol - ColPredict$Col))

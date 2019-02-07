# The following script explores the results by Gardner et al. 2017 genome 
# reserach that were downloaded from ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/
load('D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData')
load('D:/L1polymORF/Data/L1RefRanges_hg19.Rdata')

# Read in vcf file with variant calls
MEVarCall <- read.table("D:/L1polymORF/Data/nstd144.GRCh37.variant_call.vcf", 
                        as.is = T,
                        col.names = c("Chrom", "Pos", "ID", "Alt", "Type", "V6", 
                                      "V7", "Info"))
MEVarCall$Info[1:10]
MEVarRegion <- read.table("D:/L1polymORF/Data/nstd144.GRCh37.variant_region.vcf", as.is = T)
table(MEVarRegion$V5)
MEVarRegion$V8[1]
# Read sample data
SampleInfo <- read.delim("D:/L1polymORF/Data/igsr_samples.tsv", as.is = T)
table(SampleInfo$Data.collections)
SampleInfo$Data.collections[1:5]

# Get index of LINE1 deletions and check whether an allele frequency is counted
idxDel <- grep("<DEL:ME:LINE1>", MEVarCall$Type)
MEVarCall$Info[idxDel[1:10]]
grep("AF=", MEVarCall$Info[idxDel])

idxIns <- grep("<INS:ME:LINE1>", MEVarCall$Type)
grep("AF=", MEVarCall$Info[idxIns])
length(idxIns)
idxIns[1]
MEVarCall$Info[idxIns[1]]

# Extract allele frequency from info column
x <- MEVarCall$Info[idxIns[1]]
GetAF <- function(x){
  xSplit <- strsplit(x, ";")[[1]]
  AFch   <- strsplit(xSplit[length(xSplit)], "=")[[1]][2]
  as.numeric(AFch)
}
MEVarCall$Freq <- sapply(MEVarCall$Info, GetAF)

GetInsLength <- function(x){
  xSplit <- strsplit(x, ";")[[1]]
  SVLENPart <- grep("SVLEN=", xSplit, value = T)
  SVLEN   <- strsplit(SVLENPart, "=")[[1]][2]
  as.numeric(SVLEN)
}
MEVarCall$InsLength <- sapply(MEVarCall$Info, GetInsLength)
hist(MEVarCall$InsLength[idxIns])

MEVarCall$Info[which(MEVarCall$InsLength < 0)[1]]

# Create GRanges object for MEVarCall
MEVarCall$ChromName <- paste("chr", MEVarCall$Chrom, sep = "")
MEVar_GR <- makeGRangesFromDataFrame(df = MEVarCall,
                                     seqnames.field = "ChromName",
                                     start.field = "Pos",
                                     end.field = "Pos")

# Find overlap between the two GRanges and plot frequency estimates against 
# each other
OL_1000G_MEVar <- findOverlaps(L1_1000G_GR_hg19, MEVar_GR)
plot(L1_1000G_reduced$Frequency[OL_1000G_MEVar@from], 
     MEVarCall$Freq[OL_1000G_MEVar@to])
lines(c(0, 1), c(0, 1))

plot(log(L1_1000G_reduced$Frequency[OL_1000G_MEVar@from]), 
     log(MEVarCall$Freq[OL_1000G_MEVar@to]), col = rgb(0, 0, 0, alpha = 0.2),
     pch = 16)
lines(c(-10, 1), c(-10, 1), col = "red")
FreqSmu <- supsmu(log(L1_1000G_reduced$Frequency[OL_1000G_MEVar@from]), 
     log(MEVarCall$Freq[OL_1000G_MEVar@to]))
lines(FreqSmu$x, FreqSmu$y, col = "blue")

NewGR <- resize(MEVar_GR, 100, fix = "center")
sum(overlapsAny(L1GRanges,  NewGR))
OL_Ref_MEVar <- findOverlaps(L1GRanges,  MEVar_GR)
L1GRanges[OL_Ref_MEVar@from[1:10]]
MEVar_GR[OL_Ref_MEVar@to[1:10]]

#############################################################
#                                                           #
#   Plot frequency quantiles vs. insertion length           #
#         1000 genome data                                  #
#                                                           #
#############################################################

# Create a vector of L1 start classes
L1_1000G$InsLengthClass <- cut(L1_1000G$InsLength, breaks = 
                                 seq(0, 6500, 500))

# Get mean L1 frequency per start
L1WidthAggregated <- aggregate(L1_1000G[,c("InsLength", "Frequency")], 
                               by = list(L1_1000G$InsLengthClass), FUN = mean)
L1WidthAggregated_var <- aggregate(L1_1000G[,c("InsLength", "Frequency")], 
                                   by = list(L1_1000G$InsLengthClass), FUN = var)
L1WidthAggregated_n <- aggregate(L1_1000G[,c("InsLength", "Frequency")], 
                                 by = list(L1_1000G$InsLengthClass), FUN = length)

# Create matrix of quantiles per length mat
LengthClasses <- unique(L1_1000G$InsLengthClass)

ProbV <- seq(0, 1, 0.1)
QMat <- sapply(LengthClasses, function(x){
  idxL1 <- which(L1_1000G$InsLengthClass == x)
  quantile(L1_1000G$Frequency[idxL1], ProbV)
})
QMat <- QMat[,-ncol(QMat)]
QMatLog <- log(QMat)
QMatLog[is.infinite(QMatLog)] <- -10
dim(L1WidthAggregated)
Cols <- rainbow(nrow(QMatLog))
par( mfrow = c(1, 1), oma = c( 0.2,  0.2,  0.2,  0.2), 
     mai = c(1, 1, 0.2, 1),
     cex.lab = 1)
i <- 1
plot(L1WidthAggregated$InsLength, QMatLog[i,], ylim = c(-10, 0),
     type = "l", col = Cols[i])
for (i in 2:nrow(QMatLog)){
  lines(L1WidthAggregated$InsLength, QMatLog[i,], col = Cols[i])
}

CreateDisplayPdf('D:/L1polymORF/Figures/FreqVsL1Width_quantiles.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 5, width = 5)


#############################################################
#                                                           #
#   Plot frequency quantiles vs. insertion length           #
#         1000 genome data                                  #
#                                                           #
#############################################################

# Create a vector of L1 start classes
MEVarCall$InsLengthClass[idxIns] <- cut(MEVarCall$InsLength[idxIns], breaks = 
                                 seq(0, 6500, 500))

# Get mean L1 frequency per start
L1WidthAggregated <- aggregate(MEVarCall[idxIns,c("InsLength", "Freq")], 
                               by = list(MEVarCall$InsLengthClass[idxIns]), FUN = mean)
L1WidthAggregated_var <- aggregate(MEVarCall[idxIns,c("InsLength", "Freq")], 
                                   by = list(MEVarCall$InsLengthClass[idxIns]), FUN = var)
L1WidthAggregated_n <- aggregate(MEVarCall[idxIns,c("InsLength", "Freq")], 
                                 by = list(MEVarCall$InsLengthClass[idxIns]), FUN = length)

# Create matrix of quantiles per length mat
LengthClasses <- unique(MEVarCall$InsLengthClass[idxIns])

ProbV <- seq(0, 1, 0.25)
QMat <- sapply(LengthClasses, function(x){
  idxL1 <- which(MEVarCall$InsLengthClass == x)
  quantile(MEVarCall$Freq[idxL1], ProbV, na.rm = T)
})
as.numeric(MEVarCall$Freq)
QMat <- QMat[,-ncol(QMat)]
QMatLog <- log(QMat)
QMatLog[is.infinite(QMatLog)] <- -10
dim(L1WidthAggregated)
Cols <- rainbow(nrow(QMatLog))
par( mfrow = c(1, 1), oma = c( 0.2,  0.2,  0.2,  0.2), 
     mai = c(1, 1, 0.2, 1),
     cex.lab = 1)
i <- 1
plot(L1WidthAggregated$InsLength, QMatLog[i,], ylim = c(-10, 0),
     type = "l", col = Cols[i])
for (i in 2:nrow(QMatLog)){
  lines(L1WidthAggregated$InsLength, QMatLog[i,], col = Cols[i])
}

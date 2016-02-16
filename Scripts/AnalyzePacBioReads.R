# Load results
list.files("D:/L1polymORF/Data/")
load("D:/L1polymORF/Data/ReadsPerL1.RData")

ScannedReads[[1]]$pos
NrMapped <- sapply(ScannedReads, function(x) sum(!is.na(x$pos)))
hist(NrMapped, breaks= seq(0, 20000, 10),xlim = c(0, 1000))
which(NrMapped == 0)


dim(CoverMat)
CoverInL1  <- rowSums(CoverMat[,7000:13000])
CoverOutL1 <- rowSums(CoverMat[,-c(7000:13000)])
CoverInL1[which(NrMapped == 0)]
CoverOutL1/CoverInL1
hist(CoverInL1, breaks = seq(0, 290000, 1000))
min(CoverInL1)
sum(CoverInL1 < )

plot(NrMapped, CoverInL1, xlim = c(0, 1000))
plot(NrMapped, CoverOutL1/CoverInL1, xlim = c(0, 1000),
     ylim = c(0, 1))

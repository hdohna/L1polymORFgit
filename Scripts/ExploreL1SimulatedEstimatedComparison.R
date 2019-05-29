load("D:/L1polymORF/Data/L1Simulated_MELT.RData")

# Add numeric columns for L1 start and end
L1_1000G$L1StartNum <- as.numeric(as.character(L1_1000G$L1Start))
L1_1000G$L1EndNum   <- as.numeric(as.character(L1_1000G$L1End))

plot(L1_1000G$InsLength, L1_1000G$L1EndNum  - L1_1000G$L1StartNum)
lines(c(0, 6000), c(0, 6000), col = "red", lwd = 2)
plot(L1Detect$L1widthTrue, L1Detect$L1widthEst)
plot(L1Detect_Group$L1widthTrue, L1Detect_Group$L1widthEst)
hist(L1Detect$L1widthTrue)
hist(L1Detect$L1widthEst)
table(L1Detect$L1widthEst)
CurrentChrom = "chr1"
ChrT <- CreateInsertTxt(1:5)

L1Detect[which(L1Detect$SampleID == "HG00107")[1:10], ]


L1_1000G[1:5 ,c("CHROM", "POS", "L1Start", "L1End")]


par(mfrow = c(1, 2), oma = c(3, 4, 1, 1),
    mar = )
plot(1:10, 1:10, ylab = "", xlab = "X label")
plot(1:10, 1:10, ylab = "", xlab = "X label")
mtext(side = 2, "Y label", outer = T)


load('D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData')
load('D:/L1polymORF/Data/L1RefRanges_hg19.Rdata')

findOverlaps(L1_1000G_GR_hg19, L1GRanges)
L1_1000G_GR_hg19[2529]
L1GRanges[736]
L1_1000G_GR_hg19_Plus200 <- resize(L1_1000G_GR_hg19, 500, fix = "center")
L1_1000G_GR_hg19_Plus200
findOverlaps(L1_1000G_GR_hg19_Plus200, L1GRanges)
L1_1000G_GR_hg19[1098]
L1GRanges[67]

NrHitPerWidth <- sapply(seq(100, 10000, 100), function(x){
  NewGR <- resize(L1_1000G_GR_hg19, x, fix = "center")
  sum(overlapsAny(NewGR, L1GRanges))
})
plot(NrHitPerWidth)

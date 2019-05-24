library(seqinr)
L1Consens <- read.fasta("D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fas")
MatchMat <- sapply(L1Consens$L1HS_L1_Homo_sapiens, function(x){
  x == L1Consens$L1HS_L1_Homo_sapiens
})
image(MatchMat[1:100, 1:100])
MatchMat[1:10, 1:10]
dotPlot(L1Consens$L1HS_L1_Homo_sapiens[1:100], L1Consens$L1HS_L1_Homo_sapiens[1:100], 
        wsize = 3, wstep = 3)
dotPlot(L1Consens$L1HS_L1_Homo_sapiens[c(1:500, 5500:6064)], 
        L1Consens$L1HS_L1_Homo_sapiens[c(1:500, 5500:6064)])
length(L1Consens$L1HS_L1_Homo_sapiens)

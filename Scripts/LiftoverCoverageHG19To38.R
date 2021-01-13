# Load necessary packages
library(biglm)
library(Rsamtools)
library(rtracklayer)

# Load data on L1 coverage
load("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1CoverageResults.RData")

# Create genomic ranges from table that has coverage of L1 positions
L1CoverTable$Chromosome <- paste("chr", L1CoverTable$Chromosome, sep = "")
L1Cover_GR <- makeGRangesFromDataFrame(L1CoverTable, 
                                       seqnames.field = "Chromosome",
                                       start.field = "Pos",
                                       end.field = "Pos")
L1Cover_GR_HG38_List <- liftOver(L1Cover_GR, 
                                 chain = import.chain("D:/OneDrive - American University of Beirut/L1polymORF/Data/hg19ToHg38.over.chain"))
NMatched <- sapply(L1Cover_GR_HG38_List, length)
L1Cover_GR_HG38 <- unlist(L1Cover_GR_HG38_List)
L1CoverTable <- L1CoverTable[NMatched == 1, ]

# Load genome
library(BSgenome.Hsapiens.UCSC.hg19)

# Get table with L1 repeats
RepeatTable <- read.delim("D:/L1polymORF/Data/repeatL1hg19")

chr1Seq <- BSgenome.Hsapiens.UCSC.hg19[['chr1']]

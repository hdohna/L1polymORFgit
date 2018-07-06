# Process chain file
ChainLinesWithComments <- readLines("D:/L1polymORF/Data/hg19.hg19.all.chain")
FirstChar <- substr(ChainLinesWithComments, 1, 1)
sum(FirstChar == "#")
writeLines(ChainLinesWithComments[FirstChar != "#"],
           "D:/L1polymORF/Data/hg19.hg19.processed.chain")

# The following script processes the S1 table by Belyaeva et al. 2017 (PNAS)

# Load packages
library(GenomicRanges)

# Read in the data table
BelyaevaTable <- read.csv("D:/L1polymORF/Data/Belyaeva2017TableS1.csv",
                          as.is = T)

# Turn cluster into genomic ranges
ClusterStarts <- unlist(lapply(BelyaevaTable$Clusters, function(x) 
  strsplit(x, ",")[[1]]))
ClusterTable <- t(sapply(ClusterStarts, function(x) strsplit(x, ":")[[1]]))
ClusterTable <- as.data.frame(ClusterTable)
colnames(ClusterTable) <- c("chrom", "start")
ClusterTable$start <- as.numeric(ClusterTable$start)
ClusterTable$end <- ClusterTable$start + 250000
BelyaevaClustersGR <- makeGRangesFromDataFrame(ClusterTable)

# Save genomic ranges
save(list = "BelyaevaClustersGR", file = "D:/L1polymORF/Data/BelyaevaClustersGR.RData")


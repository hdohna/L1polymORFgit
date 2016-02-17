# The following script reads a table from the dbRIP database
dbRIP <- read.delim("D:/L1polymORF/Data/L1_hg19_v2h.txt", header = F)
dbRIP <- dbRIP[,-24]
colnames(dbRIP)<- c("bin", "chrom", "chromStart", "chromEnd", "name", "score",
                    "strand","originalId", "forwardPrimer", "reversePrimer", 
                    "polyClass", "polyFamily", "polySubfamily", "polySeq", 
                    "polySource", "reference", "ascertainingMethod", "remarks",
                    "tm", "fillsize", "emptysize", "disease", "genoRegion")
table(dbRIP[,12])

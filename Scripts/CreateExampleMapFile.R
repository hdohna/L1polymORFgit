# The following script identifies columns without proper genotype notation


# 
MEI1000GLines <- readLines("/labs/dflev/hzudohna/1000Genomes/ALL.chr10.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf", 
                           n = 300)
StartChar  <- substr(MEI1000GLines, 1, 1)
StartLines <- MEI1000GLines[StartChar == "#"]

ColNames  <- MEI1000GLines[max(which(StartChar == "#"))]
ColNames  <- gsub("#", "", ColNames)
ColNames  <- strsplit(ColNames, "\t")[[1]]
colnames(L1Table) <- ColNames

# Read in L1 table
L1Table <- read.delim("/labs/dflev/hzudohna/1000Genomes/LINE1all.vcf", sep = "\t",
                      skip = 2)
# Get indices of columns that contain 2
blnCol2 <- apply(L1Table, 2, FUN = function(x) length(grep("2", as.character(x))) > 0)

# Create a map table
MapTable <- cbind(L1Table[1:5, c("CHROM", "ID")], GenPos = 0.01, L1Table[1:5, c("POS")])

# Write out an example map file
write.table(MapTable, file = "/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/Chr10ExampleMap",
              row.names = F, col.names = F, quote = F, sep = " ")

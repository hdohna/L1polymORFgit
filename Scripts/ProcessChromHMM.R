# The following script processes the files with regulatory elements for 
# different cell types, assembled by the Broad Institute 
# Reference: Ernst J, Kheradpour P, Mikkelsen TS, Shoresh N, Ward LD, 
#            Epstein CB, Zhang X, Wang L, Issner R, Coyne M et al. Mapping and
#            analysis of chromatin state dynamics in nine human cell types. 
#            Nature. 2011 May 5;473(7345):43-9

# Data source: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/

# Load necessary packages
library(GenomicRanges)

# Specify folder with files (one per cell type)
DataFolder <- "D:/L1polymORF/Data/EncodeBroadHMM/"
OutputFile <- "ChromHMMcombined.txt"

# Get all file names
AllFiles <- list.files(DataFolder, full.names = T)
AllFiles <- AllFiles[-grep("txt.", AllFiles)]
AllFiles <- AllFiles[-grep(OutputFile, AllFiles)]

# Read in table with regulatory elements
RegTableCombined <- read.table(AllFiles[1], 
                      col.names = c("bin", "chrom", "start",
                                    "end", "name", "score", "strand",  
                                     "thickStart", "thickEnd", "itemRgb"))
CType <- strsplit(strsplit(AllFiles[1], "Hmm")[[1]][2], "HMM")[[1]][1]
RegTableCombined$CellType <- CType

# Loop over files and add entries that do not overlap 
for (i in 2:length(AllFiles)){
  cat("Processing file", i, "of", length(AllFiles), "\n")
  # Read in new table
  RegTableNew <- read.table(AllFiles[i],
                   col.names = c("bin", "chrom", "start", "end", "name", 
                                 "score","strand", "thickStart", "thickEnd", 
                                 "itemRgb"))

  # Determine cell type of new file
  CType <- strsplit(strsplit(AllFiles[i], "Hmm")[[1]][2], "HMM")[[1]][1]
  RegTableNew$CellType <- CType
  
  # Make genomic ranges and determine overlap
  GR1 <- makeGRangesFromDataFrame(RegTableCombined)
  GR2 <- makeGRangesFromDataFrame(RegTableNew)
  OL <- findOverlaps(GR1, GR2)
  blnNameMatch <- RegTableCombined$name[OL@from] == 
    RegTableNew$name[OL@to]
  
  # Add cell type to the matching ranges
  idxBoth <- unique(OL@from[which(blnNameMatch)])
  RegTableCombined$CellType[idxBoth] <- paste(
    RegTableCombined$CellType[idxBoth], CType, sep = ",")
  
  # Get indices of entries to append
  idxNew <- setdiff(1:nrow(RegTableNew), OL@to[blnNameMatch])
  RegTableCombined <- rbind(RegTableCombined, RegTableNew[idxNew,])
  cat(length(idxNew), "new entries added \n")
}

# Write out resulting file
write.table(RegTableCombined, 
          file = paste(DataFolder, OutputFile, sep = ""),
          row.names = F)

# Display how often certain cell type combinations occur
ComboCount <- table(RegTableCombined$CellType)
sort(ComboCount)

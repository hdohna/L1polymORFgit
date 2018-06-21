# Read info data on sources for DNAse data
DNAseInfoShort <- read.delim('D:/L1polymORF/Data/EncodeDNAse/wgEncodeRegDnaseClusteredSources.txt',
                        header = F, col.names = c("ID", "cell.1"),
                        as.is = T)

# Read data on encode cell types
EncodeInfo <- read.csv('D:/L1polymORF/Data/EncodeCellTypes.csv', as.is = T)

# Merge cell type info with DNAse source info
DNAseInfo <- merge(DNAseInfoShort, EncodeInfo)

# Group cell type IDs into stem and not stem cells
idxStem1   <- grep("stem", DNAseInfo$Tissue.5)
idxStem2   <- grep("stem", DNAseInfo$Description.3)
idxCancer  <- grep("cancer", DNAseInfo$Karyotype)
StemIDs    <- as.character(DNAseInfo$ID[union(idxStem1, idxStem2)])
CancerIDs  <- as.character(DNAseInfo$ID[idxCancer])
NotStemIDs <- as.character(setdiff(DNAseInfo$ID, union(CancerIDs, StemIDs)))

# Read DNAse hypersensitivity data
DNAseData  <- read.delim('D:/L1polymORF/Data/EncodeDNAse/wgEncodeRegDnaseClusteredV3.txt', 
                header = F, 
                col.names = c("bin", "chrom", "start", "end", "name", "score",
                              "sourceCount", "sourceIds", "sourceScores"))

# Summarize score by stem-cells and non-stem cells (move to separate script)
ScoreByStem <- t(sapply(1:nrow(DNAseData), function(i){
  IDs <- strsplit(as.character(DNAseData$sourceIds[i]), ",")[[1]]
  Scores <- as.numeric(strsplit(as.character(DNAseData$sourceScores[i]), ",")[[1]])
  c(StemScores = sum(Scores[IDs %in% StemIDs]),
    NotStemScores = sum(Scores[IDs %in% NotStemIDs]))
}))

# Save data 
save.image('D:/L1polymORF/Data/DNAseInfo.RData')

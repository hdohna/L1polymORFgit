# The following script compares L1 locations based on 1000 genome data, PacBio
# capture data and PacBio MtSinai data

load("D:/L1polymORF_old/Data/L1NonReference_Pacbio_Round2.Rdata")

# Read table with L1 info from capture experiment
FullL1Info <- read.csv("D:/L1polymORF/Data/L1InsertionInfo.csv")
FullL1GR <- makeGRangesFromDataFrame(FullL1Info)
# Get indices from files
idxFromBam <- sapply(FileNames, function(x){
  PathSplit <- strsplit(x, "/")[[1]]
  Fname <- PathSplit[length(PathSplit)]
  as.numeric(strsplit(Fname, "_")[[1]][2])
})

# Select ranges
L1RangesPacBio <- IslGRanges_reduced[idxFromBam]

# Read in 1000 Genome L1 insertion (created in script Create_1000G_L1GRanges.R)
load("D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData")

# Calculate distance from each L1 to closest L1 in 1000 genome
Dist2_1000GL1 <- distanceToNearest(L1RangesPacBio, L1_1000G_GR_hg19_NA12878, ignore.strand = T) 
Dist2_Capture<- distanceToNearest(L1RangesPacBio, FullL1GR, ignore.strand = T) 
blnClose <- Dist2_1000GL1@elementMetadata@listData$distance <= 50
sum(Dist2_Capture@elementMetadata@listData$distance <= 50)
plot(colMeans(CoverMat), type = "s")
blnFullL <- apply(CoverMat, 1, FUN =function(x) all(x > 0))

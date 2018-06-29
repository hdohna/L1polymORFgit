# The following script reads in a a bed file of replication phases and determines 
# whether full-length L1 are overrepresented in a particular phase

library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)
library(Hotelling)

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.r')

############
#  Process L1 data
############

# Read repeat table with L1HS in hg19
RepeatTable <- read.table("D:/L1polymORF/Data/repeatsHg19_L1HS", header = T)

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
L1RefGR <- GRanges(seqnames = RepeatTable$genoName,
                   ranges = IRanges(start = RepeatTable$genoStart,
                                    end = RepeatTable$genoEnd),
                   strand = RepeatTable$strand)

# Path to L1 catalogue file 
L1Catalogue <- read.csv("D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv", as.is = T)


L1Catalogue$Allele[is.na(L1Catalogue$Allele)] <- 1

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
blnInRef1         <- (L1Catalogue$end_HG38 - L1Catalogue$start_HG38) > 6000 
L1CatalogL1Mapped <- L1Catalogue[blnL1Mapped & blnAllele1,]

# Create genomic ranges for catalog L1
L1CatalogGR <- GRanges(seqnames = L1CatalogL1Mapped$Chromosome,
                       ranges = IRanges(start = pmin(L1CatalogL1Mapped$start_HG38,
                                                     L1CatalogL1Mapped$end_HG38),
                                        end = pmax(L1CatalogL1Mapped$start_HG38,
                                                   L1CatalogL1Mapped$end_HG38)),
                       strand = L1CatalogL1Mapped$strand_L1toRef)
L1CatGR_hg19 <- liftOver(L1CatalogGR, 
                         chain = import.chain("D:/L1polymORF/Data/hg38ToHg19.over.chain"))
NrRanges <- sapply(L1CatGR_hg19, length)
idxUniqueMapped <- NrRanges == 1
L1CatGR_hg19 <- unlist(L1CatGR_hg19[idxUniqueMapped])
L1CatGR_hg19 <- L1CatGR_hg19[width(L1CatGR_hg19) > 6000]

# Find overlap with fragments and catalog elements
L1RefGR_fragm    <- L1RefGR[width(L1RefGR) < 6000]

# Get combined genomic ranges and create a bed file
L1GR_combined <- c(L1CatGR_hg19, L1RefGR_fragm)
L1GR_combinedDF <- as.data.frame(L1GR_combined)
write.table(L1GR_combinedDF[1:1000,1:3], "D:/L1polymORF/Data/L1GRtable1",
            quote = F, row.names = F, col.names = F)
write.table(L1GR_combinedDF[1001:nrow(L1GR_combinedDF),1:3], "D:/L1polymORF/Data/L1GRtable2",
            quote = F, row.names = F, col.names = F)
write.table(L1GR_combinedDF[,1:3], "D:/L1polymORF/Data/L1GRtableFull",
            quote = F, row.names = F, col.names = F)

############
#  Process histone data
############

# Read table with open compartment data
HistoTableFiles <- list.files("D:/L1polymORF/Data/", 
   pattern = "H3K", full.names = T)

# Function to make a histone table
FileName <- HistoTableFiles[1]
makeHistoTable <- function(FileName, L1GR_full = L1CatGR_hg19, 
                           L1GR_fragm = L1RefGR_fragm) {
  
  # Read lines and create table with histone values and genomicRanges
  HistoLines <- readLines(FileName)
  NameSplt   <- strsplit(FileName, "/")[[1]]
  HistName   <- NameSplt[length(NameSplt)]
  HistName   <- strsplit(HistName, "per")[[1]][1]
  blnChr     <- substr(HistoLines, 1, 3) == "chr"
  HistoTable <- read.delim(text = HistoLines[blnChr], header = F,
                           col.names = c("chromosome", "start", "end", HistName))
  HistoGR    <- makeGRangesFromDataFrame(HistoTable)

  # Aggregate histone values per full-length L1   
  OverlapsFull  <- findOverlaps(HistoGR, L1GR_full)
  HistAgg_Full  <- aggregate(
    HistoTable[OverlapsFull@from, HistName] ~ OverlapsFull@to, FUN = mean)
  colnames(HistAgg_Full) <- c("idxGR", HistName)
  HistAgg_Full$start <- start(L1GR_full)[HistAgg_Full$idxGR]
  HistAgg_Full$end   <- end(L1GR_full)[HistAgg_Full$idxGR]
  HistAgg_Full$chromosome <- as.vector(seqnames(L1GR_full))[HistAgg_Full$idxGR]
  HistAgg_Full$L1State    <- "full"
  
  # Aggregate histone values per fragment L1   
  OverlapsFragm  <- findOverlaps(HistoGR, L1RefGR_fragm)
  HistAgg_Fragm  <- aggregate(
    HistoTable[OverlapsFragm@from,HistName] ~ OverlapsFragm@to, FUN = mean)
  colnames(HistAgg_Fragm)  <- c("idxGR", HistName)
  HistAgg_Fragm$start      <- start(L1GR_fragm)[HistAgg_Fragm$idxGR]
  HistAgg_Fragm$end        <- end(L1RefGR_fragm)[HistAgg_Fragm$idxGR]
  HistAgg_Fragm$chromosome <- as.vector(seqnames(L1RefGR_fragm))[HistAgg_Fragm$idxGR]
  HistAgg_Fragm$L1State    <- "fragm"
  
  # Put info for full-length and fragment L1
  rbind(HistAgg_Full, HistAgg_Fragm)
}

# Create an aggregated replication table
HistAgg <- makeHistoTable(HistoTableFiles[1])
for (FileName in HistoTableFiles[-1]){
  NewTable <- makeHistoTable(FileName)
  HistAgg  <- merge(HistAgg, NewTable)
}

############
#  Test for differences betw full-length and fragment L1
############

# Correlation between different marks
cor(HistAgg[,c("H3K27Ac", "H3K4Me1", "H3K4Me3")])
boxplot(log(HistAgg$H3K4Me3) ~ HistAgg$L1State)

# T tests
t.test(HistAgg$H3K27Ac ~ HistAgg$L1State)
t.test(HistAgg$H3K4Me1 ~ HistAgg$L1State)
t.test(HistAgg$H3K4Me3 ~ HistAgg$L1State)

# Test based on random sampling of fragments
HistVals <- HistAgg$H3K4Me3
NrFull <- sum(HistAgg$L1State == "full")
MeanFull <- mean(HistVals[HistAgg$L1State == "full"])
Samples <- SampleQuantiles(HistVals[HistAgg$L1State == "fragm"], NrFull)
sum(Samples$SampleMeans <= MeanFull) / length(Samples$SampleMeans)
  

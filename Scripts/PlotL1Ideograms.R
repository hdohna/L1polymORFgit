# The following script plots ideograms and L1 locations on them
library(GenomicRanges)
library(IdeoViz)

# Get ideo table
hg38ideo <- getIdeo("hg38")

# Read in table with known L1 
L1Catalog <- read.csv("D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv", 
                        as.is = T)
L1Catalog$values <- 1

# Get numeric activity 
ActivityNum <- L1Catalog$Activity
ActivityNum <- gsub("<", "", ActivityNum)
ActivityNum <- as.numeric(ActivityNum)
L1Catalog$ActivityNum <- ActivityNum

# Create genomic ranges
blnL1Mapped       <- !is.na(L1Catalog$start_HG38) & L1Catalog$Allele == 1 
  
L1CatalogMapped   <- L1Catalog[blnL1Mapped,]
L1CatalogMapped$newStart_HG38   <- pmin(L1CatalogMapped$start_HG38, 
                                         L1CatalogMapped$end_HG38)
L1CatalogMapped$newEnd_HG38     <- pmax(L1CatalogMapped$start_HG38, 
                                         L1CatalogMapped$end_HG38)
L1CatalogReduced <- L1CatalogMapped[,c("Chromosome", "newStart_HG38", "newEnd_HG38", "values")]
L1GR  <- makeGRangesFromDataFrame(L1CatalogReduced, seqnames.field = "Chromosome", 
                         start.field = "newStart_HG38", end.field ="newEnd_HG38",
                         keep.extra.columns = T)
L1GR  <- GenomicRanges::resize(L1GR, 100000)
L1Chr <- unique(as.vector(seqnames(L1GR)))

# Make ideo plot
ChrSet1 <- paste("chr", 1:11, sep = "")
ChrSet2 <- paste("chr", c(12:23, "X", "Y"), sep = "")
plotOnIdeo(L1Chr[L1Chr %in% ChrSet1], ideoTable = hg38ideo, 
           values_GR = L1GR, plotType = "rect",
           lwd = 20, xaxt = "n", yaxt = "n", addScale = F, col = "red")

dev.copy2pdf(file = "D:/L1polymORF/Figures/L1Ideo1_11_red.pdf")
plotOnIdeo(L1Chr[L1Chr %in% ChrSet2], ideoTable = hg38ideo, 
           values_GR = L1GR, plotType = "rect",
           lwd = 20, xaxt = "n", yaxt = "n", addScale = F, col = "red")
dev.copy2pdf(file = "D:/L1polymORF/Figures/L1Ideo12_X_red.pdf")

# Ideoplot with black lines
plotOnIdeo(L1Chr[L1Chr %in% ChrSet1], ideoTable = hg38ideo, 
           values_GR = L1GR, plotType = "rect",
           lwd = 20, xaxt = "n", yaxt = "n", addScale = F, col = "black")

dev.copy2pdf(file = "D:/L1polymORF/Figures/L1Ideo1_11_black.pdf")
plotOnIdeo(L1Chr[L1Chr %in% ChrSet2], ideoTable = hg38ideo, 
           values_GR = L1GR, plotType = "rect",
           lwd = 20, xaxt = "n", yaxt = "n", addScale = F, col = "black")
dev.copy2pdf(file = "D:/L1polymORF/Figures/L1Ideo12_X_black.pdf")

# Get only genomic ranges of Seleme L1s
L1CatalogMapped <- subset(L1Catalog, subset = Allele == 2)
L1CatalogReduced <- L1CatalogMapped[,c("Chromosome", "start_HG38", "end_HG38", "values")]
L1GR  <- makeGRangesFromDataFrame(L1CatalogReduced, seqnames.field = "Chromosome", 
                                  start.field = "start_HG38", end.field ="end_HG38",
                                  keep.extra.columns = T)
L1GR  <- GenomicRanges::resize(L1GR, 100000)
L1Chr <- unique(as.vector(seqnames(L1GR)))

# Ideoplot with black lines
plotOnIdeo(L1Chr, ideoTable = hg38ideo, 
           values_GR = L1GR, plotType = "rect",
           lwd = 20, xaxt = "n", yaxt = "n", addScale = F, col = "black")

dev.copy2pdf(file = "D:/L1polymORF/Figures/L1Ideo_Seleme.pdf")




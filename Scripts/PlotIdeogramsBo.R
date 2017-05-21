# The following script plots ideograms and L1 locations on them
library(GenomicRanges)
library(IdeoViz)

# Get ideo table
hg19ideo <- getIdeo("hg19")

# Read in bed files with insertions
AluL1    <- read.delim("D:/L1polymORF/Data/K562.highconf.novel.Alu.L1.bed")
AluL1$values    <- 1
AluL1$Alu    <- 1*(AluL1$Type == "Alu")
AluL1$L1     <- 1*(AluL1$Type == "L1")
AluL1Reduced <- AluL1[,c("Chromosome", "Left_boundary", "Right_boundary", "values")]
AluL1GR  <- makeGRangesFromDataFrame(AluL1, start.field = "Left_boundary",
                                     end.field = "Right_boundary",keep.extra.columns = T)
AluL1GR  <- GenomicRanges::resize(AluL1GR, 100000)
AluL1Chr <- unique(as.vector(seqnames(AluL1GR)))
table(AluL1$Chromosome)
table(as.vector(seqnames(AluL1GR)))
Cols <- rainbow(4)[c(2,4)]

ColorSuffix <- "Rainbow"

# Define 2 chromosome sets
ChrSet1 <- paste("chr", 1:11, sep = "")
ChrSet2 <- paste("chr", c(12:23, "X", "Y"), sep = "")

# Ideoplot with two different colors
plotOnIdeo(AluL1Chr[AluL1Chr %in% ChrSet1], ideoTable = hg19ideo,
           values_GR = AluL1GR, value_cols = c("Alu","L1"), plotType = "rect",
           xaxt = "n", yaxt = "n", addScale = F, col = Cols, val_range = c(0, 2))
dev.copy2pdf(file = paste("D:/L1polymORF/Figures/AluL1Ideo2Colors_Bo1_", ColorSuffix,
  ".pdf", sep = ""))

plotOnIdeo(AluL1Chr[AluL1Chr %in% ChrSet2], ideoTable = hg19ideo, 
           values_GR = AluL1GR, value_cols = c("Alu", "L1"), plotType = "rect",
           xaxt = "n", yaxt = "n", addScale = F, col = Cols, val_range = c(0, 2))
dev.copy2pdf(file = paste("D:/L1polymORF/Figures/AluL1Ideo2Colors_Bo2_", ColorSuffix,
             ".pdf", sep = ""))


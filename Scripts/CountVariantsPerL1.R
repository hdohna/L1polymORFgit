# The following script counts the number of variants per L1

# Source start script
#source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load necessary packages
library(Rsamtools)
library(rtracklayer)

# Read genomic ranges of LINE-1
#L1GR <-  import.bed("/labs/dflev/hzudohna/RefSeqData/L1Ranges.bed")
L1GR <-  import.bed("D:/L1polymORF/Data/L1HSRefRanges_hg19.bed")
blnFull <- width(L1GR) >= 6000

# Read vcf with variants in LINE-1s 
#L1Variants <- ReadVCF("/labs/dflev/hzudohna/1000Genomes/VariantsInL1.vcf.recode.vcf")
L1Variants <- ReadVCF("D:/L1polymORF/Data/VariantsInL1.vcf.recode.vcf")

# Create a GRanges object
L1VariantGR <- makeGRangesFromDataFrame(L1Variants, seqnames.field = "X.CHROM",
                                      start.field = "POS",
                                      end.field = "POS")

# Count the number of variants per L1
L1VarCount <- countOverlaps(L1GR, L1VariantGR)

# Perform poisson regression 
VarCountGLM <- glm(L1VarCount ~ width(L1GR) + blnFull, 
                   family = poisson(link = "identity"))
summary(VarCountGLM)
# VarCountGLM <- glm(L1VarCount ~ width(L1GR) + blnFull, 
#                    family = poisson(link = "log"))
# summary(VarCountGLM)

L1widthOrder <- order(width(L1GR))
L1VarCountSmoothed <- supsmu(width(L1GR), L1VarCount)
plot(width(L1GR), L1VarCount, col = rgb(0, 0, 0, alpha = 0.2), xlab = "L1 length",
     ylab = "Number variants")
lines(width(L1GR)[L1widthOrder], predict(VarCountGLM)[L1widthOrder],
      col = "red")
plot(width(L1GR), L1VarCount / width(L1GR), col = rgb(0, 0, 0, alpha = 0.2), xlab = "L1 length",
     ylab = "Number variants")
boxplot(L1VarCount / width(L1GR) ~ blnFull, xlab = "Full-length L1",
     ylab = "Number variants per bp")


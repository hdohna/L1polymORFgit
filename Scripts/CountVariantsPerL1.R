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
blnPlus <- as.vector(strand(L1GR) == "+")
blnPlusFull <- blnPlus[blnFull]

# Get UTR5, ORF1, ORF2, and UTR3, of full-length L1
GR_UTR5_full <- L1GR[blnFull]
GR_UTR5_full[blnPlusFull]  <- resize(L1GR[blnFull & blnPlus], 908)
GR_UTR5_full[!blnPlusFull] <- resize(L1GR[blnFull & (!blnPlus)], 908, fix = "end")

GR_ORF1_full <- L1GR[blnFull]
GR_ORF1_full[blnPlusFull]  <- narrow(L1GR[blnFull & blnPlus], start = 909, end = -4126)
GR_ORF1_full[!blnPlusFull] <- narrow(L1GR[blnFull & (!blnPlus)], start = 4126, end = -909)

GR_ORF2_full <- L1GR[blnFull]
GR_ORF2_full[blnPlusFull]  <- narrow(L1GR[blnFull & blnPlus], start = 1988, end = -236)
GR_ORF2_full[!blnPlusFull] <- narrow(L1GR[blnFull & (!blnPlus)], start = 236, end = -1988)

GR_UTR3_full <- L1GR[blnFull]
GR_UTR3_full[blnPlusFull]  <- resize(L1GR[blnFull & blnPlus], 235, fix = "end")
GR_UTR3_full[!blnPlusFull] <- resize(L1GR[blnFull & (!blnPlus)], 235)

# Read vcf with variants in LINE-1s 
#L1Variants <- ReadVCF("/labs/dflev/hzudohna/1000Genomes/VariantsInL1.vcf.recode.vcf")
L1Variants <- ReadVCF("D:/L1polymORF/Data/VariantsInL1.vcf.recode.vcf")

# Create a GRanges object
L1VariantGR <- makeGRangesFromDataFrame(L1Variants, seqnames.field = "X.CHROM",
                                      start.field = "POS",
                                      end.field = "POS")

# Count the number of variants per L1
L1VarCount <- countOverlaps(L1GR, L1VariantGR)

L1VarCount_UTR5 <- countOverlaps(GR_UTR5_full, L1VariantGR)
L1VarCount_UTR3 <- countOverlaps(GR_UTR3_full, L1VariantGR) 
L1VarCount_ORF1 <- countOverlaps(GR_ORF1_full, L1VariantGR) 
L1VarCount_ORF2 <- countOverlaps(GR_ORF2_full, L1VariantGR) 


plot(L1VarCount[blnFull], L1VarCount_UTR5 + L1VarCount_UTR3 +
       L1VarCount_ORF1 + L1VarCount_ORF2)

L1VarCountPerRange <- rbind(
  data.frame(Region = "UTR5", Count = L1VarCount_UTR5, 
             Width = width(GR_UTR5_full),
             stringsAsFactors = F),
  data.frame(Region = "UTR3", Count = L1VarCount_UTR3, 
             Width = width(GR_UTR3_full),
             stringsAsFactors = F),
  data.frame(Region = "ORF1", Count = L1VarCount_ORF1, 
             Width = width(GR_ORF1_full),
             stringsAsFactors = F),
  data.frame(Region = "ORF2", Count = L1VarCount_ORF2, 
             Width = width(GR_ORF2_full),
             stringsAsFactors = F))
  
boxplot(Count/Width ~ Region, data = L1VarCountPerRange)
VarPerRegion <- lm(Count/Width ~ Region, data = L1VarCountPerRange)
anova(VarPerRegion)

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


# The following script counts the number of variants per L1

# Source start script
#source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load necessary packages
library(Rsamtools)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

# Set the start of ORF1, ORF2 and L1 width. the values below were obtained by
# submitting L1 consensus to L1Xplorer http://l1base.charite.de/l1xplorer.php
startORF1 <- 908
startORF2 <- 1988
start3UTR <- 5812
endL1     <- 6047

# Read repeat masker table for L1HS
L1Table <- read.csv("D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv", as.is = T)
L1Table$ChrNr <- substr(L1Table$genoName, 4, nchar(L1Table$genoName))
L1Table$genoStartMinus1000 <- L1Table$genoStart - 1000
L1Table$genoEndPlus1000    <- L1Table$genoEnd + 1000
L1Table$idx <- 1:nrow(L1Table)

# Determine genomic start and end for each feature of L1 (5'UTR, ORF1, ORF2 and
# 3' UTR) on L1 on the plus strand
L1Table$start5UTR <- L1Table$genoStart
L1Table$end5UTR   <- pmin(L1Table$genoStart + startORF1 - L1Table$repStart, 
                         L1Table$genoEnd)
L1Table$startORF1 <- L1Table$genoStart + pmax(0, startORF1 - L1Table$repStart)
L1Table$endORF1   <- pmin(L1Table$genoStart + startORF2 - L1Table$repStart, 
                         L1Table$genoEnd)
L1Table$startORF2 <- L1Table$genoStart + pmax(0, startORF2 - L1Table$repStart)
L1Table$endORF2   <- pmin(L1Table$genoStart + start3UTR - L1Table$repStart, 
                         L1Table$genoEnd)
L1Table$start3UTR <- L1Table$genoStart + pmax(0, start3UTR - L1Table$repStart)
L1Table$end3UTR   <- pmin(L1Table$genoStart + start3UTR - L1Table$repStart, 
                         L1Table$genoEnd)

# Determine genomic start and end for each feature of L1 (5'UTR, ORF1, ORF2 and
# 3' UTR) on L1 on the minus strand
blnMinus <- L1Table$strand == "-"
L1Table$end5UTR[blnMinus]   <- L1Table$genoEnd[blnMinus]
L1Table$start5UTR[blnMinus] <- (L1Table$genoEnd - startORF1 + 
                                              L1Table$repLeft)[blnMinus]
L1Table$endORF1[blnMinus]   <- pmin(L1Table$genoEnd, L1Table$genoEnd - startORF1 + 
                                L1Table$repLeft)[blnMinus]
L1Table$startORF1[blnMinus] <- pmax(L1Table$genoStart, L1Table$genoEnd - startORF2 + 
                                  L1Table$repLeft)[blnMinus]
L1Table$endORF2[blnMinus]   <- pmin(L1Table$genoEnd, L1Table$genoEnd - startORF2 + 
                                     L1Table$repLeft)[blnMinus]
L1Table$startORF2[blnMinus] <- pmax(L1Table$genoStart, L1Table$genoEnd - start3UTR + 
                                  L1Table$repLeft)[blnMinus]
L1Table$end3UTR[blnMinus]   <- pmin(L1Table$genoEnd, L1Table$genoEnd - start3UTR + 
                                     L1Table$repLeft)[blnMinus]
L1Table$start3UTR[blnMinus] <- pmax(L1Table$genoStart, L1Table$genoEnd - endL1 + 
                                  L1Table$repLeft)[blnMinus]

# Read genomic ranges of LINE-1
#L1GR <-  import.bed("/labs/dflev/hzudohna/RefSeqData/L1Ranges.bed")
# L1GR    <-  import.bed("D:/L1polymORF/Data/L1HSRefRanges_hg19.bed")

# Read vcf with variants in LINE-1s 
L1Variants  <- ReadVCF("D:/L1polymORF/Data/VariantsInL1.recode.vcf")
L1Var_Left  <- ReadVCF("D:/L1polymORF/Data/VariantsInL1_leftFlank.recode.vcf")
L1Var_Right <- ReadVCF("D:/L1polymORF/Data/VariantsInL1_rightFlank.recode.vcf")
L1Variants$chromosome <- paste("chr", L1Variants$X.CHROM, sep = "")
L1Var_Left$chromosome <- paste("chr", L1Var_Left$X.CHROM, sep = "")
L1Var_Right$chromosome <- paste("chr", L1Var_Right$X.CHROM, sep = "")

# Create GRanges objects with L1 Seqences
L1GR <- makeGRangesFromDataFrame(L1Table, seqnames.field = "genoName",
                                      start.field = "genoStart",
                                      end.field = "genoEnd")
blnFull <- width(L1GR) >= 6000
blnPlus <- as.vector(strand(L1GR) == "+")
blnPlusFull <- blnPlus[blnFull]

# Get ranges of L1s and their flanks
L1GR_left <- makeGRangesFromDataFrame(L1Table, 
                                      seqnames.field = "genoName",
                                      start.field = "genoStartMinus1000",
                                      end.field = "genoStart")

L1GR_right <- makeGRangesFromDataFrame(L1Table, 
                                               seqnames.field = "genoName",
                                               start.field = "genoEnd",
                                               end.field = "genoEndPlus1000")

# Get ranges for UTR5, ORF1, ORF2, and UTR3, of full-length L1
blnUTR5 <- L1Table$start5UTR <= L1Table$end5UTR
blnORF1 <- L1Table$startORF1 <= L1Table$endORF1
blnORF2 <- L1Table$startORF2 <= L1Table$endORF2
blnUTR3 <- L1Table$start3UTR <= L1Table$end3UTR
any(L1Table$endUTR5 > L1Table$genoEnd)
any(L1Table$endORF1 > L1Table$genoEnd)
any(L1Table$blnORF2 > L1Table$genoEnd)
any(L1Table$blnUTR3 > L1Table$genoEnd)
sum(L1Table$start5UTR < L1Table$genoStart)
sum((L1Table$start5UTR - L1Table$genoStart)[L1Table$start5UTR < L1Table$genoStart])
any(L1Table$startORF1 < L1Table$genoStart)
any(L1Table$start3UTR < L1Table$genoStart)
any(L1Table$start3UTR < L1Table$genoStart)
UTR5_GR <- makeGRangesFromDataFrame(
             L1Table[blnUTR5, c("genoName", "start5UTR", "end5UTR", "strand", "idx")], 
             seqnames.field = "genoName", start.field = "start5UTR",
             end.field = "end5UTR", keep.extra.columns = T)
ORF1_GR <- makeGRangesFromDataFrame(
  L1Table[blnORF1, c("genoName", "startORF1", "endORF1", "strand", "idx")], 
  seqnames.field = "genoName", start.field = "startORF1",
  end.field = "endORF1", keep.extra.columns = T)
ORF2_GR <- makeGRangesFromDataFrame(
  L1Table[blnORF2, c("genoName", "startORF2", "endORF2", "strand", "idx")], 
  seqnames.field = "genoName", start.field = "startORF2",
  end.field = "endORF2", keep.extra.columns = T)
UTR3_GR <- makeGRangesFromDataFrame(
  L1Table[blnUTR3, c("genoName", "start3UTR", "end3UTR", "strand", "idx")], 
  seqnames.field = "genoName", start.field = "start3UTR",
  end.field = "end3UTR", keep.extra.columns = T)

  
# Get UTR5, ORF1, ORF2, and UTR3, of full-length L1
GR_UTR5_full <- L1GR[blnFull]
GR_UTR5_full[blnPlusFull]  <- resize(L1GR[blnFull & blnPlus], 908)
GR_UTR5_full[!blnPlusFull] <- resize(L1GR[blnFull & (!blnPlus)], 908, fix = "end")

GR_ORF1_full <- L1GR[blnFull]
GR_ORF1_full[blnPlusFull]  <- narrow(L1GR[blnFull & blnPlus], start = 909, 
                                     end = -4126)
GR_ORF1_full[!blnPlusFull] <- narrow(L1GR[blnFull & (!blnPlus)], start = 4126, 
                                     end = -909)

GR_ORF2_full <- L1GR[blnFull]
GR_ORF2_full[blnPlusFull]  <- narrow(L1GR[blnFull & blnPlus], start = 1988, 
                                     end = -236)
GR_ORF2_full[!blnPlusFull] <- narrow(L1GR[blnFull & (!blnPlus)], start = 236, 
                                     end = -1988)

GR_UTR3_full <- L1GR[blnFull]
GR_UTR3_full[blnPlusFull]  <- resize(L1GR[blnFull & blnPlus], 235, fix = "end")
GR_UTR3_full[!blnPlusFull] <- resize(L1GR[blnFull & (!blnPlus)], 235)

# Create a GRanges object
L1VarGR <- makeGRangesFromDataFrame(L1Variants, 
                                      start.field = "POS",
                                      end.field = "POS")
L1VarGR_Left <- makeGRangesFromDataFrame(L1Var_Left, 
                                    start.field = "POS",
                                    end.field = "POS")
L1VarGR_Right <- makeGRangesFromDataFrame(L1Var_Right, 
                                    start.field = "POS",
                                    end.field = "POS")

# Count the number of variants per L1
L1VarCount <- countOverlaps(L1GR, L1VarGR)
L1VarCount_Left  <- countOverlaps(L1GR_left, L1VarGR_Left)
L1VarCount_Right <- countOverlaps(L1GR_right, L1VarGR_Right)
L1VarCount_Flank <- L1VarCount_Left + L1VarCount_Right

cor.test(L1VarCount_Left, L1VarCount_Right) 

# Create data frame that keeps track for each triplet whether it contains a
# SNP, the type of triplet, the L1 region, the neighboring SNP density
# and whether the L1 is full-length or fragment
Seq_UTR5 <- getSeq(BSgenome.Hsapiens.UCSC.hg19, UTR5_GR[1])
TriFreq  <- trinucleotideFrequency(Seq_UTR5[[1]]) 
TriNucRC <- sapply(names(TriFreq), function(x) {
  DNASt <- DNAString(x)
  as.character(reverseComplement(DNASt))
})
TriNucMatch <- match(names(TriFreq), TriNucRC)
idxLeft     <- 1:length(TriNucMatch)
idxUnique   <- NULL
while (length(idxLeft) > 0){
  idxUnique <- c(idxUnique, idxLeft[1])
  idxLeft <- setdiff(idxLeft, c(idxLeft[1], TriNucMatch[idxLeft[1]]))
}

# Put the data together
cat("Putting the three datasets together ...\n\n")
PutDataTogether <- function(GR, L1Region, TriNucMatch = TriNucMatch, idxUnique = idxUnique){
  SeqSet <- getSeq(BSgenome.Hsapiens.UCSC.hg19, GR)
  TriFreq  <- c()
  for(i in which(width(GR) >= 3)){
    TriFreq_local <- trinucleotideFrequency(SeqSet[[i]]) 
    TriFreq_local <- TriFreq_local + TriFreq_local[TriNucMatch]
    TriFreq       <- c(TriFreq, TriFreq_local[idxUnique])
  }
  idxGR_Unrep    <- rep(GR@elementMetadata@listData$idx[width(GR) >= 3],
                             each = length(idxUnique))
  idxGR <-   rep(idxGR_Unrep, TriFreq)
  data.frame(TriNames = rep(names(TriFreq), TriFreq),
             idxGR =    idxGR,
             blnSNP = 0,
             VarCount_Flank = L1VarCount_Flank[idxGR],
             L1Region = L1Region,
             blnFull = blnFull[idxGR])
}
SNPInfo_UTR5 <- PutDataTogether(GR = UTR5_GR, L1Region = "UTR5", 
                                TriNucMatch = TriNucMatch, idxUnique = idxUnique)
cat("Building data for ORF1 ...\n")
SNPInfo_ORF1 <- PutDataTogether(ORF1_GR, "ORF1", 
                                TriNucMatch = TriNucMatch, idxUnique = idxUnique)
cat("Building data for ORF2 ...\n")
SNPInfo_ORF2 <- PutDataTogether(ORF2_GR, "ORF2", 
                                TriNucMatch = TriNucMatch, idxUnique = idxUnique)
SNPInfo_UTR3 <- PutDataTogether(UTR3_GR, "UTR3", 
                                TriNucMatch = TriNucMatch, idxUnique = idxUnique)
SNPInfo <- rbind(SNPInfo_UTR5, SNPInfo_ORF1, SNPInfo_ORF2, SNPInfo_UTR3)
cat("... done!\n")

# Paste combination of index and trinuclotide names to 
TriNuc_idx0 <- paste(SNPInfo$TriNames, SNPInfo$idxGR)

# Resize variants to get trinucleotides
L1VarGR_Tri        <- resize(L1VarGR, 3, fix = "center")
L1Var_TriNucSeq    <- getSeq(BSgenome.Hsapiens.UCSC.hg19, L1VarGR_Tri)
L1Var_TriNucSeq_RC <- reverseComplement(L1Var_TriNucSeq)
L1Var_TriNuc       <- as.character(L1Var_TriNucSeq)
blnNoMatch <- ! L1Var_TriNuc %in% names(TriFreq)[idxUnique]
sum(blnNoMatch)
length(L1Var_TriNuc)
L1Var_TriNuc[blnNoMatch] <- as.character(L1Var_TriNucSeq_RC)[blnNoMatch]
OL_L1Var      <- findOverlaps(L1GR, L1VarGR)
TriNuc_idx1   <- paste(L1Var_TriNuc[OL_L1Var@to], OL_L1Var@from)
TriMatch      <- match(TriNuc_idx1, TriNuc_idx0)
TriNuc_idx1[is.na(TriMatch)][1:100]
which(is.na(TriMatch))[1:10]

sum(is.na(TriMatch)) / length(TriMatch)
# Replace SNP indicator by one for all SNPs
TriNuc2replace <- TriNuc_idx1[!is.na(TriMatch)]
idx2ReplLeft   <- 1:length(TriNuc_idx0)
idx2Replace    <- c()
Counter <- 0
while (any(duplicated(TriNuc2replace)) & Counter < 100){
  cat(sum(duplicated(TriNuc2replace)), "duplicated entries\n")
  idx2ReplaceLocal <- unique(match(TriNuc2replace, TriNuc_idx0))
  idx2ReplaceRev   <- match(TriNuc_idx0[idx2ReplaceLocal], TriNuc2replace)
  TriNuc2replace   <- TriNuc2replace[-unique(idx2ReplaceRev)]
  idx2Replace <- c(idx2Replace, idx2ReplaceLocal)
  Counter <- Counter + 1
}
length(idx2Replace)
sum(!is.na(TriMatch))
SNPInfo$blnSNP[idx2Replace] <- 1

# Perform analysis
# SNPLogReg <- glm(blnSNP ~  TriNames + VarCount_Flank + L1Region + blnFull,
#     data = SNPInfo, family = binomial)


# Count overlaps per region
L1VarCount_UTR5 <- countOverlaps(GR_UTR5_full, L1VarGR)
L1VarCount_UTR3 <- countOverlaps(GR_UTR3_full, L1VarGR) 
L1VarCount_ORF1 <- countOverlaps(GR_ORF1_full, L1VarGR) 
L1VarCount_ORF2 <- countOverlaps(GR_ORF2_full, L1VarGR) 


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
  
par(mfrow = c(2, 2))
boxplot(L1VarCount / width(L1GR) ~ blnFull, names = c("fragment", "full-length"),
        ylab = "Variants/bp", main = "Comparison")

boxplot(Count/Width ~ Region, data = L1VarCountPerRange, ylab = "Variants/bp",
        main = "Full-length L1")
CreateDisplayPdf('D:/L1polymORF/Figures/VariantCounts.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7, width = 7)

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


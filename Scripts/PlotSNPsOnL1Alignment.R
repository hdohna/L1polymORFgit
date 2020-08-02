# The following file reads an alignment of L1 sequences on hg19 and maps hg19 
# SNPs onto that alignment

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Source required packages
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(seqinr)

# Sequence of start and end of ORF1 and ORF2
StartSeqORF1 <- "ATGGGGAAA"
StartSeqORF2 <- "ATGACAGGA"
EndSeqORF1   <- "GCCAAAATGTAA"
EndSeqORF2   <- "GGTGGGAATTGA"

# Load data on L1 coverage
load("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1CoverageResults.RData")

# Create genomic ranges
L1CoverTable$Chromosome <- paste("chr", L1CoverTable$Chromosome, sep = "")
L1Cover_GR <- makeGRangesFromDataFrame(L1CoverTable, 
                                       seqnames.field = "Chromosome",
                                       start.field = "Pos",
                                       end.field = "Pos")
# Read repeat masker table for L1HS
L1Table <- read.csv("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HS_repeat_table_Hg19.csv", as.is = T)

# Create GRanges objects with L1 Seqences
L1GR <- makeGRangesFromDataFrame(L1Table, seqnames.field = "genoName",
                                 start.field = "genoStart",
                                 end.field = "genoEnd")

# Get sequences and create a character of sequence names
L1Seq    <- getSeq(BSgenome.Hsapiens.UCSC.hg19, L1GR)
SeqNames <- paste(as.vector(seqnames(L1GR)), start(L1GR), end(L1GR), sep = "_")

# Form different subsets and write them out as fasta files
L1Aligned <- read.fasta(file = "file:///D:/OneDrive - American University of Beirut/L1polymORF/Data/L1seqHg19_minLength6000_aligned.txt")
# L1Aligned <- read.alignment(file = "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1seqHg19_minLength6000_aligned.txt",
#                             format = "fasta")

# Read vcf with variants in LINE-1s 
L1Variants  <- ReadVCF("D:/OneDrive - American University of Beirut/L1polymORF/Data/VariantsInL1.recode.vcf")
L1Variants$chromosome <- paste("chr", L1Variants$X.CHROM, sep = "")

# Create a GRanges object of variants inside L1s and their flanking regions
L1VarGR <- makeGRangesFromDataFrame(L1Variants, 
                                    start.field = "POS",
                                    end.field = "POS")

# Get indicies of sequenced positions
idxSeqPosList <- lapply(L1Aligned, function(x) which(x != "-"))

# Identify start and end of ORF1 and ORF2 on alignment
i <- 1
ORFStartEndList <- lapply(1:length(L1Aligned), function(i){
  Seq <- L1Aligned[[i]]
  idxSeq <-  which(Seq != "-")
  ORF1StartMatch <- matchPattern(StartSeqORF1, paste(toupper(Seq[idxSeq]), collapse = ""))
  ORF1EndMatch   <- matchPattern(EndSeqORF1, paste(toupper(Seq[idxSeq]), collapse = ""))
  ORF2StartMatch <- matchPattern(StartSeqORF2, paste(toupper(Seq[idxSeq]), collapse = ""))
  ORF2EndMatch   <- matchPattern(EndSeqORF2, paste(toupper(Seq[idxSeq]), collapse = ""))
  list(ORF1Start = idxSeq[start(ORF1StartMatch)],
       ORF1End = idxSeq[end(ORF1EndMatch)],
       ORF2Start = idxSeq[start(ORF2StartMatch)],
       ORF2End = idxSeq[end(ORF2EndMatch)])
  
})
table(unlist(lapply(ORFStartEndList, function(x) x$ORF1Start)))
table(unlist(lapply(ORFStartEndList, function(x) x$ORF1End)))
table(unlist(lapply(ORFStartEndList, function(x) x$ORF2Start)))
table(unlist(lapply(ORFStartEndList, function(x) x$ORF2End)))
ORF1Start <- 1682
ORF1End   <- 2733
ORF2Start <- 2798
ORF2End   <- 6990

# Get genomic ranges of L1 in the alignment
L1GRNames <-  paste(as.vector(seqnames(L1GR)), start(L1GR), end(L1GR), sep = "_")
L1Match  <- match(names(L1Aligned), L1GRNames)
L1GRFull <- L1GR[L1Match]

# Determine the start of the L1s depending on strand 
idxMinus <- 1 + c(as.vector(strand(L1GRFull)) == "-")
StartEnd <- cbind(start(L1GRFull), end(L1GRFull))
L1Starts <- sapply(seq_along(idxMinus), function(i) StartEnd[i, idxMinus[i]])

# Function to find a set of genomic positions (specified as genomic ranges 
# PosGR) on an alignment
GenPos2AlignPos <- function(PosGR){
  
  # Overlap between two sets of genomic ranges
  OL_Pos   <- findOverlaps(L1GRFull, PosGR)
  
  # Determine the positions relative to the L1 start
  RelPos <- abs(start(L1GRFull)[OL_Pos@from] - start(PosGR)[OL_Pos@to]) + 1
  
  # Determine the positions on the alignment
  idxPosDf <- data.frame() 
  for(x in unique(OL_Pos@from)){
    Seq    <- L1Aligned[[x]]
    idxSeq <-  which(Seq != "-")
    blnL1  <- OL_Pos@from == x
    PosSeq <- RelPos[blnL1]
    NewDf <- data.frame(PosInAlign = idxSeq[PosSeq],
                        idxInGR    = OL_Pos@to[blnL1])
    idxPosDf <- rbind(idxPosDf, NewDf)
  }
  idxPosDf
  # unlist(lapply(unique(OL_Pos@from), function(x){
  #   Seq    <- L1Aligned[[x]]
  #   idxSeq <-  which(Seq != "-")
  #   PosSeq <- RelPos[OL_Pos@from == x]
  #   idxSeq[PosSeq]
  # }))
}

# Get positions of SNPs in alignment
SNPposDF <- GenPos2AlignPos(L1VarGR)
SNPposAlign <- SNPposDF$PosInAlign

# Get positions of coverage data in alignment
blnCoverInL1  <- overlapsAny(L1Cover_GR, L1GRFull)
L1CoverSubset <- L1CoverTable[blnCoverInL1,]
CoverPosDF    <- GenPos2AlignPos(L1Cover_GR[blnCoverInL1])

L1CoverSubset$CoverSt <- (L1CoverSubset$CoverMean - min(L1CoverSubset$CoverMean))/
  (max(L1CoverSubset$CoverMean) - min(L1CoverSubset$CoverMean))

# Average coverage per alignment position
CoverMeanPerAlign <- aggregate(L1CoverSubset$CoverMean[CoverPosDF$idxInGR], 
                              by = list(CoverPosDF$PosInAlign),
                              FUN = mean)
CoverPosCount <- table(CoverPosDF$PosInAlign)
CountMatch <- match(CoverMeanPerAlign$Group.1, names(CoverPosCount))
CoverMeanPerAlign$Count <- CoverPosCount[CountMatch]

# Determine proportion of blanks per position
L1IndelMat <- t(sapply(L1Aligned, function(x) x == "-"))
PropIndel <- colMeans(L1IndelMat)


# PLOT SNPs and TFB Regions on L1
layout(matrix(c(1, 1, 2, 3, 4, 5), 3, 2, byrow = TRUE))
# par(oma = c(0.1,  0.2,  0.1,  0.2), 
#     mai = c(0.7, 1, 0.2, 1), cex.lab = 1.2)

par(mfrow = c(1, 1))
plot(c(-300, max(SNPposAlign)), c(0, 1), type= "n", xlab = "", ylab="",xaxt="n",frame=F,
     yaxt = "n")
segments(1, 0.05, max(SNPposAlign), 0.05) # UTRs
rect(c(ORF1Start, ORF2Start), c(0, 0), c(ORF1End, ORF2End), 0.1, 
     border = "black", col ="lightgrey") # ORFs
text(0.5 * c(ORF1Start + ORF1End, ORF2Start + ORF2End), 0.05, c("ORF1", "ORF2"), cex = 0.75)

segments(SNPposAlign, 0.15, SNPposAlign, 0.2, col= rgb(0, 0, 0, 0.03), lwd = 0.1)
text(c(-300, -300), c(0.175, 0.275), c("SNP", "Indel"), cex = 0.75)
segments(1:ncol(L1IndelMat), 0.25, 1:ncol(L1IndelMat), 0.3, col= rgb(0, 0, 0, 0.1*(1 - PropIndel)), 
         lwd = 0.1)
lines(CoverMeanPerAlign$Group.1[CoverMeanPerAlign$Count >= 10], 
      0.1 + CoverMeanPerAlign$x[CoverMeanPerAlign$Count >= 10] / 10)
min(CoverMeanPerAlign$x[CoverMeanPerAlign$Count >= 10])
max(CoverMeanPerAlign$x[CoverMeanPerAlign$Count >= 10])
axis(2, at = 0.1 + seq(0.4, 0.8, 0.2), labels = seq(4, 8, 2))
mtext("Coverage", side = 2, line = 2,  at = 0.7)


par(page = F, mai = c(0.5, 1, 0.2, 0.1))
# boxplot(L1VarCount / width(L1GR) ~ blnFull, names = c("fragment", "full-length"),
#         ylab = "SNPsper LINE-1 bp", main = "B", xlab = "LINE-1 type")
# boxplot(Count/Width ~ Region, data = L1VarCountPerRange, ylab = "",
#         xlab = "LINE-1 region",
#         main = "C")

# Get mean number of SNPs per full-length and fragment L1
PlotMeanSNP(AggPerL1Pos_FullFrag, NameV = c("Fragment", "Full-length"),
            Main = "b", YLim = c(0, 0.025), YLab = "SNPs per LINE-1 bp", Border = NA,
            PlotP = "P < 0.0001")
PlotMeanSNP(AggPerL1Pos_ORFvsUTR, NameV = c("UTR", "ORF"),
            Main = "c", YLim = c(0, 0.025), Border = NA, PlotP = "P < 0.0001")
PlotMeanSNP(AggPerORFPos[!AggPerORFPos$blnFull, ], 
            NameV = c("synonymous", "non-synonymous"),
            YLim = c(0, 0.02), Main = "d", YLab = "SNPs per LINE-1 bp", Border = NA)
PlotMeanSNP(AggPerORFPos[AggPerORFPos$blnFull, ], 
            NameV = c("synonymous", "non-synonymous"),
            YLim = c(0, 0.002), Main = "e", Border = NA,
            PlotP = "P = 0.003")

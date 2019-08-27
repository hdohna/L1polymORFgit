# The following script counts the number of variants per L1

# Source start script
#source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load necessary packages
library(biglm)
library(Rsamtools)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Load data with proportion mismatch
load("D:/L1polymORF/Data/L1HS_PropMismatch.RData")

# Load data on L1 coverage
load("D:/L1polymORF/Data/L1CoverageResults.RData")
# load("D:/FantomData/Data/hg19.cage_peak_phase1and2combined_tpm.osc_GRanges.RData")

# Create genomic ranges
L1CoverTable$Chromosome <- paste("chr", L1CoverTable$Chromosome, sep = "")
L1Cover_GR <- makeGRangesFromDataFrame(L1CoverTable, 
                                       seqnames.field = "Chromosome",
                                       start.field = "Pos",
                                       end.field = "Pos")

# Getting trinucleotide sequence around each bp within L1
L1Cover_TriGR <- resize(L1Cover_GR, width = 3, fix = "center")
cat("Getting trinucletide sequence of each bp within L1 .... ")
TriNuc        <- getSeq(BSgenome.Hsapiens.UCSC.hg19, L1Cover_TriGR)
cat("done!\n")
TriNucChar          <- as.character(TriNuc)
TriNuc_RC           <- reverseComplement(TriNuc)
TriNuc_RCChar       <- as.character(TriNuc_RC)
TriNucChar          <- pmin(TriNucChar, TriNuc_RCChar)
L1CoverTable$TriNuc <- TriNucChar

# Set the start of ORF1, ORF2 and L1 width. The values below were obtained by
# submitting L1 consensus to L1Xplorer http://l1base.charite.de/l1xplorer.php
startORF1 <- 908
startORF2 <- 1988
start3UTR <- 5812
endL1     <- 6047

# Get the start of non-synonymous
StartsNonSynORF1 <- seq(0, 1011, 3)
StartsNonSynORF2 <- seq(0, 3822, 3)

StartSeqORF1 <- "ATGGGGAAA"
StartSeqORF2 <- "ATGACAGGA"
EndSeqORF1   <- "GCCAAAATGTAA"
EndSeqORF2   <- "GGTGGGAATTGA"

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

# Boolean vector whether end of ORF2 is included
blnORF2EndIncluded <- L1Table$genoEnd > L1Table$genoStart + start3UTR - 
  L1Table$repStart + 10

# Determine start of ORF1 and 
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

# Boolean vector whether most of the 5' UTR is present
L1Table$bln5UTRPresent <- L1Table$repStart <= 300
L1Table$bln5UTRPresent[blnMinus] <- (L1Table$repLeft <= 300)[blnMinus]

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
# FantomGR_1000 <- resize(FantomGR, 1000, fix = "center")
# blnOLFantom <- overlapsAny(L1GR, FantomGR_1000, ignore.strand = T)
# sum(blnOLFantom)
ChrV <- as.vector(seqnames(L1GR))
Diff <- NULL
for(chr in unique(ChrV)){
  GRsubset <- L1GR[ChrV == chr]
  StartOrder <- order(start(GRsubset))
  NewDiff <- start(GRsubset)[StartOrder[-1]] - 
    end(GRsubset)[StartOrder[-length(StartOrder)]]
  Diff <- c(Diff, NewDiff)
}
hist(Diff, breaks = seq(-100, 3.1*10^7, 10^2),
     xlim = c(0, 10^4))
min(Diff)
sum(Diff <= 100)
L1Width <- width(L1GR)
blnFull <- width(L1GR) >= 6000
blnPlus <- as.vector(strand(L1GR) == "+")
blnPlusFull <- blnPlus[blnFull]
L1FullStart <- start(L1GR[blnFull])
L1FullEnd   <- end(L1GR[blnFull])
min(L1Width)

# Get start and end positions of ORF1 and ORF2
L1Seq            <- getSeq(BSgenome.Hsapiens.UCSC.hg19, L1GR)
ORF1StartList    <- vmatchPattern(StartSeqORF1, L1Seq, max.mismatch = 1)
ORF1EndList      <- vmatchPattern(EndSeqORF1, L1Seq, max.mismatch = 2)
blnORF1Start     <- sapply(ORF1StartList, function(x) length(x) > 0)
blnORF1End       <- sapply(ORF1EndList, function(x) length(x) > 0)
blnORF1StartEnd  <- blnORF1Start & blnORF1End
ORF2StartList    <- vmatchPattern(StartSeqORF2, L1Seq, max.mismatch = 1)
ORF2EndList      <- vmatchPattern(EndSeqORF2, L1Seq, max.mismatch = 2)
blnORF2Start     <- sapply(ORF2StartList, function(x) length(x) > 0)
blnORF2End       <- sapply(ORF2EndList, function(x) length(x) > 0)
blnORF2StartEnd  <- blnORF2Start & blnORF2End
cat(sum(blnORF1Start), "L1 with ORF1 start motif\n")
cat(sum(blnORF1End), "L1 with ORF1 end motif\n")
cat(sum(blnORF1StartEnd), "L1 with ORF1 start and end motif\n")
cat(sum(blnORF1StartEnd & blnFull), "full-length L1 with ORF1 start and end motif\n")
cat(sum(blnORF2Start), "L1 with ORF2 start motif\n")
cat(sum(blnORF2End), "L1 with ORF2 end motif\n")
cat(sum(blnORF2StartEnd), "L1 with ORF2 start and end motif\n")
cat(sum(blnORF2StartEnd & blnFull), "full-length L1 with ORF2 start and end motif\n")

# Check that all full-length LINE-1 have ends of ORF1 and 2
all(blnORF1End[blnFull])
all(blnORF2End[blnFull])
all(blnORF1Start[blnFull])
all(blnORF2Start[blnFull])

# Get starts and ends of ORF1 and ORF2 within L1
ORF1Starts    <- sapply(ORF1StartList@ends[blnORF1Start], function(x) x[1] - 8)
ORF2Starts    <- sapply(ORF2StartList@ends[blnORF2Start], function(x) x[1] - 8)
ORf1Ends      <- sapply(ORF1EndList@ends[blnORF1End], function(x) x[length(x)])
ORf2Ends      <- sapply(ORF2EndList@ends[blnORF2End], function(x) x[length(x)])
ORf1StartEnds <- sapply(which(blnORF1StartEnd), function(i) {
                   x <- ORF1StartList@ends[[i]]
                   y <- ORF1EndList@ends[[i]]
                        c(x[1] - 8, y[length(y)])
                 })
ORf2StartEnds <- sapply(which(blnORF2StartEnd), function(i) {
                  x <- ORF2StartList@ends[[i]]
                  y <- ORF2EndList@ends[[i]]
                       c(x[1] - 8, y[length(y)])
                 })
plot(ORf1Ends[which(blnORF1End) %in% which(blnFull)])
plot(ORF1Starts[which(blnORF1Start) %in% which(blnFull)])
plot(ORF1Starts[blnStartEndFull], ORf1Ends[blnStartEndFull])
plot(ORf2Ends[which(blnORF2End) %in% which(blnFull)])
plot(ORf1StartEnds[2,] - ORf1StartEnds[1,], ylim = c(1000, 1030))
table(ORf1StartEnds[2,] - ORf1StartEnds[1,])
table(ORf2StartEnds[2,] - ORf2StartEnds[1,])

# Create genomic ranges for non-synonymous and synonymous coding positions
ORF1Length    <- ORf1StartEnds[2,] - ORf1StartEnds[1,]
ORF2Length    <- ORf2StartEnds[2,] - ORf2StartEnds[1,]
idxProperORF1 <- which(blnORF1StartEnd)[ORF1Length %in% c(1013, 1016, 1022)]
idxProperORF2 <- which(blnORF2StartEnd)[ORF2Length %in% c(3821, 3824, 3827, 3830)]
idxORFEnd     <- which((blnORF1End | blnORF2End))
idxORF1End    <- which(blnORF1End)
idxORF2End    <- which(blnORF2End)
idxORF1Start  <- which(blnORF1End)
idxORF2End    <- which(blnORF2End)
blnBothORFs   <- blnORF1End & blnORF2End
# blnORF2End[blnBothORFs]  <- sapply(which(blnBothORFs), function(x){
#   j <- which(idxORF1End == x)
#   k <- which(idxORF2End == x)
#   ORf2Ends[k] > ORf1Ends[j]
# })
# idxORF2End   <- which(blnORF2End)
StartNonSyn  <- NULL
StartSyn     <- NULL
ChrVCode     <- NULL
ChrV         <- as.vector(seqnames(L1GR))
blnORFFull   <- NULL
blnORFProper <- NULL
ORFType      <- NULL
L1Start      <- start(L1GR)
L1End        <- end(L1GR)
cat("Getting genomic ranges of synonymous and nonsynonymous coding positions ...")
for (i in idxORFEnd){
#for (i in union(idxProperORF1, idxProperORF2)){
  j <- which(idxORF1End == i)
  k <- which(idxORF2End == i)
  if (blnPlus[i]){ # L1 on positive strand
    ORF1Pos <- ORf1Ends[j] - StartsNonSynORF1 - 3
    ORF2Pos <- ORf2Ends[k] - StartsNonSynORF2 - 3
    blnNeg1 <- ORF1Pos < 0
    blnNeg2 <- ORF2Pos < 0
    ORF1Pos <- ORF1Pos[!blnNeg1]
    ORF2Pos <- ORF2Pos[!blnNeg2]
    NewStartNonSyn <- L1Start[i] + c(ORF1Pos, ORF2Pos)
    StartNonSyn    <- c(StartNonSyn, NewStartNonSyn)
    StartSyn       <- c(StartSyn, NewStartNonSyn + 2)
  } else { # L1 on negative strand
    ORF1Pos <- ORf1Ends[j] - StartsNonSynORF1 - 2
    ORF2Pos <- ORf2Ends[k] - StartsNonSynORF2 - 2
    blnNeg1 <- ORF1Pos < 0
    blnNeg2 <- ORF2Pos < 0
    ORF1Pos <- ORF1Pos[!blnNeg1]
    ORF2Pos <- ORF2Pos[!blnNeg2]
    NewStartNonSyn <- L1End[i] - c(ORF1Pos, ORF2Pos)
    StartNonSyn    <- c(StartNonSyn, NewStartNonSyn)
    StartSyn       <- c(StartSyn, NewStartNonSyn - 1)
  }
  
  # Update info whether ORF is complete (blnORFFull), type of ORF and chromosome
  blnORFFull  <- c(blnORFFull, rep(all(!blnNeg1), length(ORF1Pos)),
                  rep(all(!blnNeg2), length(ORF2Pos)))
  blnORFProper <- c(blnORFProper, rep(i %in% idxProperORF1, length(ORF1Pos)),
                   rep(i %in% idxProperORF2, length(ORF2Pos)))
  ORFType  <- c(ORFType, rep("ORF1", length(ORF1Pos)),
                rep("ORF2", length(ORF2Pos)))
  ChrVCode <- c(ChrVCode, rep(ChrV[i], length(NewStartNonSyn)))
}

# Create genomic ranges
GRNonSyn <- GRanges(seqnames = ChrVCode, IRanges(start = StartNonSyn,
                                                 end = StartNonSyn + 1))
GRSyn    <- GRanges(seqnames = ChrVCode, IRanges(start = StartSyn,
                                                 end = StartSyn))

cat("done!\n")

# Check that the first three letters are AA, T
getSeq(BSgenome.Hsapiens.UCSC.hg19, GRNonSyn[1])
getSeq(BSgenome.Hsapiens.UCSC.hg19, GRSyn[1])
EndSeqORF1 <- "GCCAAAATGTAA"
EndSeqORF2 <- "GGTGGGAATTGA"

getSeq(BSgenome.Hsapiens.UCSC.hg19, GRanges(seqnames = "chr1",
                                            IRanges(10^4 + 1, 10^4 + 4), strand = "+"))
getSeq(BSgenome.Hsapiens.UCSC.hg19, GRanges(seqnames = "chr1",
                                            IRanges(10^4 + 1, 10^4 + 4), strand = "-"))
# Below is the old way to get coding sequences 
L1FullSeq     <- getSeq(BSgenome.Hsapiens.UCSC.hg19, L1GR[blnFull])
ORF1StartList <- vmatchPattern(StartSeqORF1, L1FullSeq, max.mismatch = 1)
ORF1Starts    <- sapply(ORF1StartList@ends, function(x) x[1] - 8)
ORF2StartList <- vmatchPattern(StartSeqORF2, L1FullSeq, max.mismatch = 1)
ORF2Starts    <- sapply(ORF2StartList@ends, function(x) x[1] - 8)

# Create genomic ranges for non-synonymous and synonymous coding positions
# idxFull     <- which(blnFull)
# StartNonSyn <- NULL
# StartSyn    <- NULL
# ChrVCode    <- NULL
# ChrVFull    <- as.vector(seqnames(L1GR)[idxFull])
# for (i in 1:length(idxFull)){
#   if (blnPlusFull[i]){
#     NewStartNonSyn <- L1FullStart[i] + c(ORF1Starts[i] - 1 + StartsNonSynORF1,
#                         ORF2Starts[i] - 1 + StartsNonSynORF2)
#     StartNonSyn <- c(StartNonSyn, NewStartNonSyn)
#     StartSyn <- c(StartSyn, NewStartNonSyn + 2)
#   } else {
#     NewStartNonSyn <- L1FullEnd[i] - c(ORF1Starts[i] + StartsNonSynORF1,
#                          ORF2Starts[i] + StartsNonSynORF2)
#     StartNonSyn <- c(StartNonSyn, NewStartNonSyn)
#     StartSyn    <- c(StartSyn, NewStartNonSyn - 1)
#   }
#   ChrVCode <- c(ChrVCode, rep(ChrVFull[i], length(NewStartNonSyn)))
# }
# GRNonSyn_Full <- GRanges(seqnames = ChrVCode, IRanges(start = StartNonSyn,
#                                                  end = StartNonSyn + 1))
# GRSyn_Full     <- GRanges(seqnames = ChrVCode, IRanges(start = StartSyn,
#                                                  end = StartSyn))
# 
# GRNonSyn <- c(GRNonSyn_Fragm, GRNonSyn_Full)
# GRSyn    <- c(GRSyn_Fragm, GRSyn_Full)


Exons <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
Genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
Cds <- cds(TxDb.Hsapiens.UCSC.hg19.knownGene)
Promoters <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream = 5*10^3)
CdsTx <- cdsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "tx")
idxOLCdsL1 <- which(overlapsAny(CdsTx, L1GR))
NonSynStart   <- NULL
SynStart      <- NULL
ChrVCode_Gene <- NULL

# Loop over coding sequences to get coding regions and nonsynonymous positions
# in genes
for(j in idxOLCdsL1){
  x <- CdsTx[[j]]
  if (as.vector(strand(x))[1] == "-"){
    OffSet <- 1
    for (i in 1:length(x)){
      Chr <- as.vector(seqnames(x))[1]
      NonSynStartNew <- end(x[i]) - seq(OffSet, width(x[i]), 3)
      NonSynStart    <- c(NonSynStart, NonSynStartNew)
      SynStart       <- c(SynStart, NonSynStartNew - 1)
      OffSet         <- OffSet + width(x[i]) %% 3
      ChrVCode_Gene <- c(ChrVCode_Gene, rep(Chr, length(NonSynStartNew)))
    }
  }  else {
    OffSet <- 0
    for (i in 1:length(x)){
      Chr <- as.vector(seqnames(x))[1]
      NonSynStartNew <- start(x[i]) + seq(OffSet, width(x[i]), 3)
      NonSynStart    <- c(NonSynStart, NonSynStartNew)
      SynStart       <- c(SynStart, NonSynStartNew + 2)
      OffSet         <- OffSet + width(x[i]) %% 3
      ChrVCode_Gene <- c(ChrVCode_Gene, rep(Chr, length(NonSynStartNew)))
    }
  }
}
length(ChrVCode_Gene)
length(NonSynStart)

# Create genomic ranges
GRNonSyn_Gene <- GRanges(seqnames = ChrVCode_Gene, IRanges(start = NonSynStart,
                                                 end = NonSynStart + 1))
GRSyn_Gene <- GRanges(seqnames = ChrVCode_Gene, IRanges(start = SynStart,
                                                           end = SynStart))

# Check that the first three letters are AT, G
getSeq(BSgenome.Hsapiens.UCSC.hg19, GRNonSyn[1])
getSeq(BSgenome.Hsapiens.UCSC.hg19, GRSyn[1])

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

# Create a GRanges object of variants inside L1s and their flanking regions
L1VarGR <- makeGRangesFromDataFrame(L1Variants, 
                                    start.field = "POS",
                                    end.field = "POS")

# Count the number of variants per L1
L1VarCount <- countOverlaps(L1GR, L1VarGR)
L1VarCount_Left  <- countOverlaps(L1GR_left, L1VarGR_Left)
L1VarCount_Right <- countOverlaps(L1GR_right, L1VarGR_Right)
L1VarCount_Flank <- L1VarCount_Left + L1VarCount_Right
cor.test(L1VarCount_Left, L1VarCount_Right) 

# Overlaps between individual bp and L1
OL_bpL1 <- findOverlaps(L1Cover_GR, L1GR)

# Add columns with different info
L1CoverTable$blnSNP  <- overlapsAny(L1Cover_GR, L1VarGR)
L1CoverTable$blnUTR5 <- overlapsAny(L1Cover_GR, UTR5_GR)
L1CoverTable$blnORF1 <- overlapsAny(L1Cover_GR, ORF1_GR)
L1CoverTable$blnORF2 <- overlapsAny(L1Cover_GR, ORF2_GR)
L1CoverTable$blnUTR3 <- overlapsAny(L1Cover_GR, UTR3_GR)
L1CoverTable$Exons   <- overlapsAny(L1Cover_GR, Exons)
L1CoverTable$Genes   <- overlapsAny(L1Cover_GR, Genes)
L1CoverTable$Promoters   <- overlapsAny(L1Cover_GR, Promoters)
L1CoverTable$L1VarCount_Flank <- NA
L1CoverTable$L1VarCount_Flank[OL_bpL1@from] <- L1VarCount_Flank[OL_bpL1@to]
L1CoverTable$PropMismatch <- NA
L1CoverTable$PropMismatch[OL_bpL1@from] <- PropMismatch[OL_bpL1@to]
L1CoverTable$L1Width <- NA
L1CoverTable$L1Width[OL_bpL1@from] <- L1Width[OL_bpL1@to]
L1CoverTable$blnFull <- NA
L1CoverTable$blnFull[OL_bpL1@from] <- blnFull[OL_bpL1@to]
# L1CoverTable$blnOLFantom <- NA
# L1CoverTable$blnOLFantom[OL_bpL1@from] <- blnOLFantom[OL_bpL1@to]
#L1CoverTable$CodeType <- "NonCode" 
L1CoverTable$Coding <- overlapsAny(L1Cover_GR, c(GRSyn, GRNonSyn))
L1CoverTable$NonSyn <- overlapsAny(L1Cover_GR, GRNonSyn) 
L1CoverTable$NonSyn_Full <- overlapsAny(L1Cover_GR, GRNonSyn[blnORFFull])
L1CoverTable$Coding_Full <- overlapsAny(L1Cover_GR, c(GRSyn[blnORFFull], 
                                                      GRNonSyn[blnORFFull]))
L1CoverTable$NonSyn_Proper <- overlapsAny(L1Cover_GR, GRNonSyn[blnORFProper])
L1CoverTable$Coding_Proper <- overlapsAny(L1Cover_GR, c(GRSyn[blnORFProper], 
                                                      GRNonSyn[blnORFProper]))
L1CoverTable$NonSyn_Gene <- overlapsAny(L1Cover_GR, GRNonSyn_Gene)
L1CoverTable$Coding_Gene <- overlapsAny(L1Cover_GR, c(GRSyn_Gene, GRNonSyn_Gene))

L1CoverTable$bln5UTRPresent <- NA
L1CoverTable$bln5UTRPresent[OL_bpL1@from] <- L1Table$bln5UTRPresent[OL_bpL1@to]

sum(blnORFFull) / sum(blnORFProper)
# L1CoverTable$NonSyn <- L1CoverTable$NonSyn & (!L1CoverTable$NonSyn_Gene) 
# L1CoverTable$Coding <- L1CoverTable$Coding & (!L1CoverTable$Coding_Gene) 

# Perform analysis without interaction
# cat("Performing regression analysis without interaction ... ")
# SNPLogReg <- bigglm(blnSNP ~  TriNuc + L1VarCount_Flank + blnUTR5 + blnORF1 + 
#                       blnORF2 + blnFull + CoverMean +
#                       L1Width + PropMismatch + NonSyn + Coding,
#                     data = L1CoverTable, family = binomial(), chunksize = 3*10^4,
#                     maxit = 20)
# SNPLogReg_Summary <- summary(SNPLogReg)
# SNPLogReg_SumDF   <- as.data.frame(SNPLogReg_Summary$mat)
# SNPLogReg_SumDF$ExpCoef <- exp(SNPLogReg_SumDF$Coef) 
# cat("done!\n")
sum(L1CoverTable$Coding_Full)
sum(L1CoverTable$NonSyn_Full)
sum(blnORFFull)
# Perform analysis with interaction
cat("Performing regression analysis with interaction ... ")
SNPLogRegInt <- bigglm(blnSNP ~  TriNuc + L1VarCount_Flank + CoverMean +
                      L1Width + PropMismatch + Genes + Exons + Promoters +
                        blnFull + 
                        # NonSyn + Coding + Coding_Gene + NonSyn_Gene + 
                        Coding_Full + NonSyn_Full +
                        # Coding_Proper + NonSyn_Proper +
                        # Coding*Exons +
                      #Coding*Genes + #  + Coding*Promoters +
                      #NonSyn*Exons +
                      #NonSyn*Genes + # NonSyn*Promoters + 
                      Coding_Full*blnFull, # + NonSyn_Full*blnFull,
                      # Coding_Proper*blnFull + NonSyn_Proper*blnFull,
                    data = L1CoverTable, family = binomial(), chunksize = 3*10^4,
                    maxit = 20)
cat("done!\n")

# Export data frame with regression results
SNPLogRegInt_Summary <- summary(SNPLogRegInt)
SNPLogRegInt_SumDF   <- as.data.frame(SNPLogRegInt_Summary$mat)
SNPLogRegInt_SumDF$ExpCoef <- exp(SNPLogRegInt_SumDF$Coef) 
ResultPath <- "D:/L1polymORF/Data/L1VariantRegrResults2019-07-17.csv"
cat("Writing regression results to", ResultPath, "\n")
write.csv(SNPLogRegInt_SumDF, ResultPath)

# Analyze the proportion of SNPs on non-synonymous sites 
NonsynEffect <- sapply(1:length(L1GR), function(i){
  L1Subset <- L1CoverTable[OL_bpL1@from[OL_bpL1@to == i], ]
  L1Subset <- L1Subset[L1Subset$Coding, ]
  if(any(L1Subset$blnSNP)){
    print(i)
    LReg <- glm(blnSNP ~ NonSyn, data = L1Subset)
    LRegCoeffs <- summary(LReg)$coefficients
    LRegCoeffs <- rbind(LRegCoeffs, rep(NA, 4))
    LRegCoeffs[2, c("Estimate", "Std. Error", "Pr(>|t|)")]
  } else {
    c(NA, NA, NA)
  }
})
hist(NonsynEffect[1,])
which(NonsynEffect[1,] < -0.4)
L1Table[630, ]
NonsynEffect[,630]
hist(NonsynEffect[3,], breaks = seq(0, 1, 0.01))
hist(NonsynEffect[3, NonsynEffect[1,] < 0], breaks = seq(0, 1, 0.01))
hist(NonsynEffect[3, blnFull], breaks = seq(0, 1, 0.01))
hist(NonsynEffect[1, blnFull])
plot(width(L1GR), NonsynEffect[1,])
plot(width(L1GR), NonsynEffect[1,], ylim = c(-5*10^-2, 5*10^-2))
plot(width(L1GR), NonsynEffect[3,])
lines(c(0, 10000), c(0,0), col = "red")
t.test(NonsynEffect[1,] ~ blnFull)
blnNoOutlier <- NonsynEffect[1,] < 0.05 &  NonsynEffect[1,] > -0.05 
LM <- lm(NonsynEffect[1,blnNoOutlier] ~ width(L1GR)[blnNoOutlier], 
         weights = 1/NonsynEffect[2,blnNoOutlier])
summary(LM)



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



par(mfrow = c(1, 2),  mai = c(1, 0.5, 0.2, 0.1), 
    oma = c( 0.2,  2,  0.2,  0.2))
boxplot(L1VarCount / width(L1GR) ~ blnFull, names = c("fragment", "full-length"),
        ylab = "SNPs/bp", main = "A")

boxplot(Count/Width ~ Region, data = L1VarCountPerRange, ylab = "",
        main = "B")
mtext(side =2, "SNPs/bp", outer = T)
CreateDisplayPdf('D:/L1polymORF/Figures/VariantCounts.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 4, width = 7)


# # Get UTR5, ORF1, ORF2, and UTR3, of full-length L1
# GR_UTR5_full <- L1GR[blnFull]
# GR_UTR5_full[blnPlusFull]  <- resize(L1GR[blnFull & blnPlus], startORF1)
# GR_UTR5_full[!blnPlusFull] <- resize(L1GR[blnFull & (!blnPlus)], startORF1, fix = "end")
# 
# GR_ORF1_full <- L1GR[blnFull]
# GR_ORF1_full[blnPlusFull]  <- narrow(L1GR[blnFull & blnPlus], start = 909, 
#                                      end = -4126)
# GR_ORF1_full[!blnPlusFull] <- narrow(L1GR[blnFull & (!blnPlus)], start = 4126, 
#                                      end = -909)
# 
# GR_ORF2_full <- L1GR[blnFull]
# GR_ORF2_full[blnPlusFull]  <- narrow(L1GR[blnFull & blnPlus], start = startORF2, 
#                                      end = -236)
# GR_ORF2_full[!blnPlusFull] <- narrow(L1GR[blnFull & (!blnPlus)], start = 236, 
#                                      end = -startORF2)
# 
# GR_UTR3_full <- L1GR[blnFull]
# GR_UTR3_full[blnPlusFull]  <- resize(L1GR[blnFull & blnPlus], 235, fix = "end")
# GR_UTR3_full[!blnPlusFull] <- resize(L1GR[blnFull & (!blnPlus)], 235)

# L1VarGR_Left <- makeGRangesFromDataFrame(L1Var_Left, 
#                                     start.field = "POS",
#                                     end.field = "POS")
# L1VarGR_Right <- makeGRangesFromDataFrame(L1Var_Right, 
#                                     start.field = "POS",
#                                     end.field = "POS")


# Create data frame that keeps track for each triplet whether it contains a
# SNP, the type of triplet, the L1 region, the neighboring SNP density
# and whether the L1 is full-length or fragment
# Seq_UTR5 <- getSeq(BSgenome.Hsapiens.UCSC.hg19, UTR5_GR[1])
# TriFreq  <- trinucleotideFrequency(Seq_UTR5[[1]]) 
# TriNucRC <- sapply(names(TriFreq), function(x) {
#   DNASt <- DNAString(x)
#   as.character(reverseComplement(DNASt))
# })
# TriNucMatch <- match(names(TriFreq), TriNucRC)
# idxLeft     <- 1:length(TriNucMatch)
# idxUnique   <- NULL
# while (length(idxLeft) > 0){
#   idxUnique <- c(idxUnique, idxLeft[1])
#   idxLeft <- setdiff(idxLeft, c(idxLeft[1], TriNucMatch[idxLeft[1]]))
# }
# 
# # Resize variants to get trinucleotides
# L1VarGR_Tri        <- resize(L1VarGR, 3, fix = "center")
# L1Var_TriNucSeq    <- getSeq(BSgenome.Hsapiens.UCSC.hg19, L1VarGR_Tri)
# L1Var_TriNucSeq_RC <- reverseComplement(L1Var_TriNucSeq)
# L1Var_TriNuc       <- as.character(L1Var_TriNucSeq)
# blnNoMatch         <- ! L1Var_TriNuc %in% names(TriFreq)[idxUnique]
# L1Var_TriNuc[blnNoMatch] <- as.character(L1Var_TriNucSeq_RC)[blnNoMatch]
# 
# # Function to put the data together
# PutDataTogether <- function(GR, L1Region, TriNucMatch = TriNucMatch, 
#                             idxUnique = idxUnique, L1VarGR = L1VarGR,
#                             L1Var_TriNuc = L1Var_TriNuc,
#                             MaxCounter = 100){
#   
#   # Get ssequences of each genomic region
#   SeqSet <- getSeq(BSgenome.Hsapiens.UCSC.hg19, GR)
#   TriFreq  <- c()
#   
#   # Count how often each trinucleotide occurs in each genomic region
#   for(i in which(width(GR) >= 3)){
#     TriFreq_local <- trinucleotideFrequency(SeqSet[[i]]) 
#     TriFreq_local <- TriFreq_local + TriFreq_local[TriNucMatch]
#     TriFreq       <- c(TriFreq, TriFreq_local[idxUnique])
#   }
#   
#   # Create one row for each tri-nucleotide
#   idxGRRegion_Unrep <- rep(which(width(GR) >= 3), each = length(idxUnique))
#   idxGRRegion       <- rep(idxGRRegion_Unrep, TriFreq)
#   idxGR_Unrep       <- rep(GR@elementMetadata@listData$idx[width(GR) >= 3],
#                              each = length(idxUnique))
#   idxGR             <- rep(idxGR_Unrep, TriFreq)
#   SNPData <- data.frame(TriNames = rep(names(TriFreq), TriFreq),
#                idxGR =    idxGR,
#                idxGRRegion = idxGRRegion,
#                blnSNP = 0,
#                VarCount_Flank = L1VarCount_Flank[idxGR],
#                L1Region = L1Region,
#                blnFull = blnFull[idxGR],
#                L1Width = L1Width[idxGR],
#                PropMismatch = PropMismatch[idxGR])
#    
#   # Create vector of trinucleotide IDs
#   TriNucIDs0 <- paste(SNPData$TriNames, SNPData$idxGRRegion)
#   
#   # Resize variants to get trinucleotides
#   OL_L1Var      <- findOverlaps(GR, L1VarGR)
#   TriNucIDs1    <- paste(L1Var_TriNuc[OL_L1Var@to], OL_L1Var@from)
#   TriMatch      <- match(TriNucIDs1, TriNucIDs0)
# 
#   sum(is.na(TriMatch)) / length(TriMatch)
#   # Replace SNP indicator by one for all SNPs
#   TriNuc2replace <- TriNucIDs1[!is.na(TriMatch)]
#   idx2ReplLeft   <- 1:length(TriNucIDs0)
#   idx2Replace    <- c()
#   Counter <- 0
#   while (any(duplicated(TriNuc2replace)) & Counter < MaxCounter){
#     cat(sum(duplicated(TriNuc2replace)), "duplicated entries\n")
#     idx2ReplaceLocal <- unique(match(TriNuc2replace, TriNucIDs0))
#     idx2ReplaceRev   <- match(TriNucIDs0[idx2ReplaceLocal], TriNuc2replace)
#     TriNuc2replace   <- TriNuc2replace[-unique(idx2ReplaceRev)]
#     idx2Replace      <- c(idx2Replace, idx2ReplaceLocal)
#     Counter          <- Counter + 1
#   }
#   length(idx2Replace)
#   sum(!is.na(TriMatch))
#   SNPData$blnSNP[idx2Replace] <- 1
#   SNPData
#   
# }
# cat("Building data for 5' UTR ...\n")
# SNPInfo_UTR5 <- PutDataTogether(GR = UTR5_GR, L1Region = "UTR5", 
#                                 TriNucMatch = TriNucMatch, idxUnique = idxUnique,
#                                 L1VarGR = L1VarGR,
#                                 L1Var_TriNuc = L1Var_TriNuc)
# cat("Building data for ORF1 ...\n")
# SNPInfo_ORF1 <- PutDataTogether(ORF1_GR, "ORF1", 
#                                 TriNucMatch = TriNucMatch, idxUnique = idxUnique,
#                                 L1VarGR = L1VarGR,
#                                 L1Var_TriNuc = L1Var_TriNuc)
# cat("Building data for ORF2 ...\n")
# SNPInfo_ORF2 <- PutDataTogether(ORF2_GR, "ORF2", 
#                                 TriNucMatch = TriNucMatch, idxUnique = idxUnique,
#                                 L1VarGR = L1VarGR,
#                                 L1Var_TriNuc = L1Var_TriNuc)
# cat("Building data for 3' UTR ...\n")
# SNPInfo_UTR3 <- PutDataTogether(UTR3_GR, "UTR3", 
#                                 TriNucMatch = TriNucMatch, idxUnique = idxUnique,
#                                 L1VarGR = L1VarGR,
#                                 L1Var_TriNuc = L1Var_TriNuc)
# SNPInfo <- rbind(SNPInfo_UTR5, SNPInfo_ORF1, SNPInfo_ORF2, SNPInfo_UTR3)
# cat("... done!\n")
# SNPLogRegInt <- bigglm(blnSNP ~  TriNames + VarCount_Flank + L1Region + blnFull +
#                          PropMismatch + L1Region*blnFull,
#                        data = SNPInfo, family = binomial(), chunksize = 3*10^4,
#                        maxit = 20)
# 



# VarPerRegion <- lm(Count/Width ~ Region, data = L1VarCountPerRange)
# anova(VarPerRegion)
# 
# # Perform poisson regression 
# VarCountGLM <- glm(L1VarCount ~ width(L1GR) + blnFull, 
#                    family = poisson(link = "identity"))
# summary(VarCountGLM)
# VarCountGLM <- glm(L1VarCount ~ width(L1GR) + blnFull, 
#                    family = poisson(link = "log"))
# summary(VarCountGLM)

# L1widthOrder <- order(width(L1GR))
# L1VarCountSmoothed <- supsmu(width(L1GR), L1VarCount)
# plot(width(L1GR), L1VarCount, col = rgb(0, 0, 0, alpha = 0.2), xlab = "L1 length",
#      ylab = "Number variants")
# lines(width(L1GR)[L1widthOrder], predict(VarCountGLM)[L1widthOrder],
#       col = "red")
# plot(width(L1GR), L1VarCount / width(L1GR), col = rgb(0, 0, 0, alpha = 0.2), xlab = "L1 length",
#      ylab = "Number variants")
# 
save.image("D:/L1polymORF/Data/L1VariantCount.RData")
# 

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
library(seqinr)

# File paths
ResultPathCombined  <- "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1VariantRegrResults_combined_2021-01-02.csv"
ResultPathAll  <- "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1VariantRegrResults_all_2021-01-02.csv"
ResultPathFull <- "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1VariantRegrResults_full_2021-01-02.csv"
ResultPathAllPacBio  <- "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1VariantRegrResults_all_PacBio_2021-01-02.csv"
ResultPathFullPacBio <- "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1VariantRegrResults_full_PacBio_2020-11-15.csv"


######################################
#                                    #
#        Load and process data       #
#                                    #
######################################

# Load genomic ranges of TF binding sits on L1 
load("D:/OneDrive - American University of Beirut/L1Evolution/Data/L1_TFBGR.RData")

# Load genomic ranges of SNPS that differ from HWE
load("D:/OneDrive - American University of Beirut/L1polymORF/Data/SNPsDifferFromHWE.RData")

# Load data with proportion mismatch
load("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HS_PropMismatch.RData")

# Load data on L1 coverage
load("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1CoverageResults.RData")

# Read in vcf file with MELT deletion calls
MEDelCall <- ReadVCF("D:/OneDrive - American University of Beirut/L1polymORF/Data/DEL.final_comp.vcf")
MEDelCall$chromosome <- paste("chr", MEDelCall$X.CHROM, sep = "")
MEDel_GR  <- makeGRangesFromDataFrame(df = MEDelCall,
                                      start.field = "POS",
                                      end.field = "POS")

# function to get numeric genotype
GetNumericGenotype <- function(x){
  Split1 <- strsplit(x, ":")[[1]][1]
  Split2 <- strsplit(Split1, "/")[[1]]
  sum(as.numeric(Split2))
}

# Get numeric genotype of all reference L1 deletions
GTCols <- grep("L1Filtered", colnames(MEDelCall))
length(GTCols)
L1RefNumGen <- 2 - sapply(GTCols, function(x){
  sapply(1:nrow(MEDelCall), function(y) GetNumericGenotype(MEDelCall[y,x]))
})
max(L1RefNumGen, na.rm = T)

# Add columns for frequency and sample size
MEDelCall$Freq       <- rowSums(L1RefNumGen, na.rm = T)
MEDelCall$SampleSize <- apply(L1RefNumGen, 1, function(x) 2*sum(!is.na(x)))

# Create genomic ranges from table that has coverage of L1 positions
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
endORF1   <- 1921
startORF2 <- 1988
start3UTR <- 5813
endL1     <- 6047

# Get the start of non-synonymous
StartsNonSynORF1 <- seq(0, 1011, 3)
StartsNonSynORF1 <- seq(0, 1022, 3)
StartsNonSynORF2 <- seq(0, 3822, 3)
StartsNonSynORF2 <- seq(0, 3842, 3)

StartSeqORF1      <- "ATGGGGAAA"
StartSeqORF1_long <- "ATGGGGAAAAAACAGAACAGA"
StartSeqORF2      <- "ATGACAGGA"
StartSeqORF2_long <- "ATGACAGGATCAAATTCACACATA"
EndSeqORF1        <- "GCCAAAATGTAA"
EndSeqORF2        <- "GGTGGGAATTGA"

# Read repeat masker table for L1HS
L1Table <- read.csv("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HS_repeat_table_Hg19.csv", as.is = T)
L1Table$ChrNr <- substr(L1Table$genoName, 4, nchar(L1Table$genoName))
L1Table$genoStartMinus1000 <- L1Table$genoStart - 1000
L1Table$genoEndPlus1000    <- L1Table$genoEnd + 1000
L1Table$idx <- 1:nrow(L1Table)

# Determine genomic start and end for each feature of L1 (5'UTR, ORF1, ORF2 and
# 3' UTR) on L1 on the minus strand
blnMinus <- L1Table$strand == "-"

# Read vcf with variants in LINE-1s 
L1Variants  <- ReadVCF("D:/OneDrive - American University of Beirut/L1polymORF/Data/VariantsInL1_1KG.recode.vcf")
L1VariantsCNV  <- L1Variants[grep("VT=CNV", L1Variants$INFO),]
L1Variants$VT <- sapply(L1Variants$INFO, function(x) {
  Split1 <- strsplit(x, ";")[[1]]
  grep("VT=", Split1, value = T)
  })
L1Variants$INFO[1:50]
L1Variants$AlleleFreq <- sapply(L1Variants$INFO, function(x){
  InfoSplit <- strsplit(x, ";")[[1]]
  AF <- grep("AF=", InfoSplit, value = T)
  AF <- AF[-grep("_AF=", AF)]
  as.numeric(strsplit(AF, "AF=")[[1]][2])
})
L1Variants$chromosome <- paste("chr", L1Variants$X.CHROM, sep = "")
L1Variants$VarID <- paste(L1Variants$X.CHROM, L1Variants$POS)
L1Variants_Indel <- L1Variants[L1Variants$VT == "VT=INDEL",]
L1Variants  <- L1Variants[grep("VT=SNP", L1Variants$INFO),]
hist(L1Variants$AlleleFreq)
L1Var_Left  <- ReadVCF("D:/OneDrive - American University of Beirut/L1polymORF/Data/VariantsInL1_leftFlank.recode.vcf")
L1Var_Right <- ReadVCF("D:/OneDrive - American University of Beirut/L1polymORF/Data/VariantsInL1_rightFlank.recode.vcf")
L1Var_Left$chromosome <- paste("chr", L1Var_Left$X.CHROM, sep = "")
L1Var_Right$chromosome <- paste("chr", L1Var_Right$X.CHROM, sep = "")

# Read in variants from PacBio sequencing of HG002 (obtained from https://downloads.pacbcloud.com/public/publications/2019-HG002-CCS/smallvariants/)
L1VarPacBio <- ReadVCF("D:/OneDrive - American University of Beirut/L1polymORF/Data/VariantsInL1_HG002_PacBio_withInfo_GATKHC.recode.vcf")
L1VarPacBio$chromosome <- paste("chr", L1VarPacBio$X.CHROM, sep = "")
L1VarPacBio$INFO[1:5]

# Subset to retain only SNPs
L1VarPacBio <- L1VarPacBio[grep("VRT=1", L1VarPacBio$INFO), ]

# Create a GRanges object of variants inside L1s and their flanking regions
L1VarGRPacBio <- makeGRangesFromDataFrame(L1VarPacBio, 
                                    start.field = "pos",
                                    end.field = "pos")
# Read in variants from dbSNP
L1VarDbSNP <- ReadVCF("D:/OneDrive - American University of Beirut/L1polymORF/Data/VariantsInL1_dbSNP.recode.vcf")
L1VarDbSNP$chromosome <- paste("chr", L1VarDbSNP$X.CHROM, sep = "")
L1VarGRDbSNP <- makeGRangesFromDataFrame(L1VarDbSNP, 
                                          start.field = "POS",
                                          end.field = "POS")

# Read in 1000 genome variants of HG002
L1VarHG002 <- read.table("D:/OneDrive - American University of Beirut/L1polymORF/Data/SNPsInHG002_all.vcf")
L1VarHG002$chromosome <- paste("chr", L1VarHG002$X.CHROM, sep = "")

# Create a GRanges object of variants inside L1s and their flanking regions
L1VarGRHG002 <- makeGRangesFromDataFrame(L1VarHG002, 
                                    start.field = "POS",
                                    end.field = "POS")

# Create GRanges objects with L1 Seqences
L1GR <- makeGRangesFromDataFrame(L1Table, seqnames.field = "genoName",
                                      start.field = "genoStart",
                                      end.field = "genoEnd")
sum(width(L1GR) >= 6000)

# Determine intervals between consecutive L1s
ChrV <- as.vector(seqnames(L1GR))
Diff <- NULL
for(chr in unique(ChrV)){
  GRsubset <- L1GR[ChrV == chr]
  StartOrder <- order(start(GRsubset))
  NewDiff <- start(GRsubset)[StartOrder[-1]] - 
    end(GRsubset)[StartOrder[-length(StartOrder)]]
  Diff <- c(Diff, NewDiff)
}
hist(Diff, breaks = seq(-100, 3.1*10^7, 10),
     xlim = c(-50, 10^2))
min(Diff)
sum(Diff <= 100)
sum(Diff <= 100)/length(Diff)
sum(Diff <= 0)/length(Diff)
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
blnORF1Start_2   <- sapply(ORF1StartList, function(x) length(x) > 1)
blnORF1End       <- sapply(ORF1EndList, function(x) length(x) > 0)
blnORF1StartEnd  <- blnORF1Start & blnORF1End
ORF2StartList    <- vmatchPattern(StartSeqORF2, L1Seq, max.mismatch = 1)
ORF2EndList      <- vmatchPattern(EndSeqORF2, L1Seq, max.mismatch = 2)
blnORF2Start     <- sapply(ORF2StartList, function(x) length(x) > 0)
blnORF2End       <- sapply(ORF2EndList, function(x) length(x) > 0)
blnORF2StartEnd  <- blnORF2Start & blnORF2End

# Return summaries of ORF starts and ends
cat(sum(blnORF1Start), "L1 with ORF1 start motif\n")
cat(sum(blnORF1Start_2), "L1 with 2 ORF1 start motifs\n")
cat(sum(blnORF1End),   "L1 with ORF1 end motif\n")
cat(sum(blnORF1StartEnd), "L1 with ORF1 start and end motif\n")
cat(sum(blnORF1StartEnd & blnFull), "full-length L1 with ORF1 start and end motif\n")
cat(sum(blnORF1Start_2 & blnFull), "full-length L1 with 2 ORF1 start motifs\n")
cat(sum(blnORF2Start),  "L1 with ORF2 start motif\n")
cat(sum(blnORF2End),    "L1 with ORF2 end motif\n")
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
                   x <- ORF1StartList@ends[[i]] - 8
                   y <- ORF1EndList@ends[[i]]
                   c(x[1], y[length(y)])
                   Xrep <- rep(x, each = length(y))
                   Yrep <- rep(y, length(x))
                   Diff <- Yrep - Xrep
                   idxMin <- which.min(abs(Diff - 1016))

                   idxX <- 0:(length(Xrep) - 1) %/% length(y) + 1
                   idxY <- 0:(length(Xrep) - 1) %% length(x) + 1
                   c(x[idxX[idxMin]], y[idxY[idxMin]])
})
ORf2StartEnds <- sapply(which(blnORF2StartEnd), function(i) {
                  x <- ORF2StartList@ends[[i]] - 8
                  y <- ORF2EndList@ends[[i]]
                  c(x[1], y[length(y)])
                       # Xrep <- rep(x, each = length(y))
                       # Yrep <- rep(y, length(x))
                       # Diff <- Yrep - Xrep
                       # idxMin <- which.min(abs(Diff - 3827))
                       # idxX <- 0:(length(Xrep) - 1) %/% length(y) + 1
                       # idxY <- 0:(length(Xrep) - 1) %% length(x) + 1
                       # c(x[idxX[idxMin]], y[idxY[idxMin]])
})
plot(ORf1Ends[which(blnORF1End) %in% which(blnFull)])
plot(ORF1Starts[which(blnORF1Start) %in% which(blnFull)])
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
idxORF1Start  <- which(blnORF1Start)
idxORF2Start  <- which(blnORF2Start)
blnBothORFs   <- blnORF1End & blnORF2End

# Get ranges of L1s and their flanks
L1GR_left <- makeGRangesFromDataFrame(L1Table, 
                                      seqnames.field = "genoName",
                                      start.field = "genoStartMinus1000",
                                      end.field = "genoStart")

L1GR_right <- makeGRangesFromDataFrame(L1Table, 
                                       seqnames.field = "genoName",
                                       start.field = "genoEnd",
                                       end.field = "genoEndPlus1000")

# Create a GRanges object of variants inside L1s and their flanking regions
L1VarGR <- makeGRangesFromDataFrame(L1Variants, 
                                    start.field = "POS",
                                    end.field = "POS")
sum(overlapsAny(L1VarGRDbSNP, L1VarGR))

# Create a GRanges object of indels inside L1s and their flanking regions
L1VarIndelGR <- makeGRangesFromDataFrame(L1Variants_Indel, 
                                         start.field = "POS",
                                         end.field = "POS")

#####################################################
#                                                   #
#    Identify non-synonymous positions on ORFs      #
#                                                   #
#####################################################

# Initialize objects needed for loop that defines non-synonymous positions
StartNonSynInt <- NULL
StartSynInt    <- NULL
StrandV        <- NULL
ChrVCodeInt    <- NULL
ChrV           <- as.vector(seqnames(L1GR))
blnORFFull     <- NULL
blnORFProper   <- NULL
blnORFProperInt <- NULL
ORFType        <- NULL
ORFTypeInt     <- NULL
ORFPosInt      <- NULL
ORFPos0Int     <- NULL
L1Start        <- start(L1GR)
L1End          <- end(L1GR)
cat("Getting genomic ranges of synonymous and nonsynonymous coding positions ...")
i <- 1

# Loop over all L1 that have at least one ORF end motif
sum(blnPlus)/length(blnPlus)
for (i in idxORFEnd){
  j <- which(idxORF1End == i)
  k <- which(idxORF2End == i)
  l <- which(idxORF1Start == i)
  m <- which(idxORF2Start == i)
  if (blnPlus[i]){ # L1 on positive strand
    
    # Count codon starts from the start
    ORF1Pos_st <- ORF1Starts[l] + StartsNonSynORF1
    ORF2Pos_st <- ORF2Starts[m] + StartsNonSynORF2
    
    # Count codon starts from the end 
    ORF1Pos <- ORf1Ends[j] - StartsNonSynORF1 - 2
    ORF2Pos <- ORf2Ends[k] - StartsNonSynORF2 - 2
    blnNeg1 <- ORF1Pos < 0
    blnNeg2 <- ORF2Pos < 0
    ORF1Pos <- ORF1Pos[!blnNeg1]
    ORF2Pos <- ORF2Pos[!blnNeg2]
    ORF1PosInt <- intersect(ORF1Pos_st, ORF1Pos)
    ORF2PosInt <- intersect(ORF2Pos_st, ORF2Pos)
    ORF1PosIntMatch <- match(ORF1PosInt, ORF1Pos_st)
    ORF2PosIntMatch <- match(ORF2PosInt, ORF2Pos_st)
    ORF1Pos0Int     <- StartsNonSynORF1[ORF1PosIntMatch]
    ORF2Pos0Int     <- StartsNonSynORF2[ORF2PosIntMatch]
    NewStartNonSynInt <- L1Start[i] + c(ORF1PosInt, ORF2PosInt) - 1
  } else { # L1 on negative strand
    ORF1Pos_st <- ORF1Starts[l] + StartsNonSynORF1 + 1
    ORF2Pos_st <- ORF2Starts[m] + StartsNonSynORF2 + 1
    ORF1Pos    <- ORf1Ends[j] - StartsNonSynORF1 - 1
    ORF2Pos    <- ORf2Ends[k] - StartsNonSynORF2 - 1
    blnNeg1    <- ORF1Pos < 0
    blnNeg2    <- ORF2Pos < 0
    ORF1Pos    <- ORF1Pos[!blnNeg1]
    ORF2Pos    <- ORF2Pos[!blnNeg2]
    ORF1PosInt <- intersect(ORF1Pos_st, ORF1Pos)
    ORF2PosInt <- intersect(ORF2Pos_st, ORF2Pos)
    ORF1PosIntMatch <- match(ORF1PosInt, ORF1Pos_st)
    ORF2PosIntMatch <- match(ORF2PosInt, ORF2Pos_st)
    ORF1Pos0Int     <- StartsNonSynORF1[ORF1PosIntMatch]
    ORF2Pos0Int     <- StartsNonSynORF2[ORF2PosIntMatch]
    NewStartNonSynInt <- L1End[i] - c(ORF1PosInt, ORF2PosInt) + 1
  }
  
  # Update vectors
  ORFPosInt  <- c(ORFPosInt,  c(ORF1PosInt, ORF2PosInt))
  ORFPos0Int <- c(ORFPos0Int, c(ORF1Pos0Int, ORF2Pos0Int))
  StartNonSynInt <- c(StartNonSynInt, NewStartNonSynInt)
  
  # Update info whether ORF is complete (blnORFFull), type of ORF and chromosome
  blnORFFull  <- c(blnORFFull, rep(all(!blnNeg1), length(ORF1Pos)),
                  rep(all(!blnNeg2), length(ORF2Pos)))
  blnORFProperInt <- c(blnORFProperInt, rep(i %in% idxProperORF1, length(ORF1PosInt)),
                      rep(i %in% idxProperORF2, length(ORF2PosInt)))
  ORFTypeInt  <- c(ORFTypeInt, rep("ORF1", length(ORF1PosInt)),
                rep("ORF2", length(ORF2PosInt)))
  ChrVCodeInt <- c(ChrVCodeInt, rep(ChrV[i], length(NewStartNonSynInt)))
  StrandV     <- c(StrandV, rep(c("-", "+")[1 + blnPlus[i]], length(NewStartNonSynInt)))
  
}
sum(blnORFFull)
length(blnORFProperInt)
length(StartNonSynInt)
sum(blnORFProperInt)
length(blnORFProperInt)

# Create genomic ranges
GRNonSynInt1 <- GRanges(seqnames = ChrVCodeInt, 
                       IRanges(start = StartNonSynInt,
                               end = StartNonSynInt),
                       strand = StrandV,
                       mcols = data.frame(ORFPos = ORFPosInt,
                                          ORFPos0 = ORFPos0Int + (StrandV == "-"),
                                          ORFType = ORFTypeInt))
GRNonSynInt2 <- GRanges(seqnames = ChrVCodeInt, 
                        IRanges(start = StartNonSynInt + 1,
                                end = StartNonSynInt + 1),
                        strand = StrandV,
                        mcols = data.frame(ORFPos = ORFPosInt,
                                           ORFPos0 = ORFPos0Int + (StrandV == "+"),
                                           ORFType = ORFTypeInt))
GRNonSynInt <- c(GRNonSynInt1, GRNonSynInt2)
GRSynInt    <- GRanges(seqnames = ChrVCodeInt, IRanges(start = StartNonSynInt + 2,
                                                    end = StartNonSynInt + 2),
                    strand = StrandV,
                    mcols = data.frame(ORFPos = ORFPosInt,
                                       ORFPos0 = ORFPos0Int + 2,
                                       ORFType = ORFTypeInt))
# Separate genomic ranges by ORF1 and 2
GRNonSynInt_ORF1 <- GRNonSynInt[GRNonSynInt@elementMetadata@listData$mcols.ORFType == "ORF1"]
GRNonSynInt_ORF2 <- GRNonSynInt[GRNonSynInt@elementMetadata@listData$mcols.ORFType == "ORF2"]
GRSynInt_ORF1    <- GRSynInt[GRSynInt@elementMetadata@listData$mcols.ORFType == "ORF1"]
GRSynInt_ORF2    <- GRSynInt[GRSynInt@elementMetadata@listData$mcols.ORFType == "ORF2"]

#################################################################
#                                                               #
#    Check identification non-synonymous positions on ORFs      #
#                                                               #
#################################################################

# Check what proportion of starting positions of ORF1 and ORF2 correspond
# to positions based on repeat masker
GR1_ORF1 <- GRNonSynInt[GRNonSynInt@elementMetadata@listData$mcols.ORFPos0 == 0 &
  GRNonSynInt@elementMetadata@listData$mcols.ORFType == "ORF1"]
OL1 <- findOverlaps(ORF1_GR, GR1_ORF1)
sapply(-15:15, function(i) sum(start(ORF1_GR)[OL1@from] == start(GR1_ORF1)[OL1@to] + i))
cbind(start(ORF1_GR)[OL1@from],start(GR1_ORF1)[OL1@to] )
start(GR1_ORF1)[OL1@to] - start(ORF1_GR)[OL1@from]

GR1_ORF2 <- GRNonSynInt[GRNonSynInt@elementMetadata@listData$mcols.ORFPos0 == 0 &
                          GRNonSynInt@elementMetadata@listData$mcols.ORFType == "ORF2"]
OL2 <- findOverlaps(ORF2_GR, GR1_ORF2)
sapply(-15:15, function(i) sum(start(ORF1_GR)[OL1@from] == start(GR1_ORF1)[OL1@to] + i))
cbind(start(ORF2_GR)[OL2@from], start(GR1_ORF2)[OL2@to])
start(GR1_ORF2)[OL2@to] - start(ORF2_GR)[OL2@from]
sum(start(ORF2_GR)[OL2@from] == start(GR1_ORF2)[OL2@to])

# Check what proportion of position is discovered by forming intersection
# or checking for appropriate ORF length
# length(GRNonSynInt) / length(GRNonSyn)
# sum(blnORFFull) / length(GRNonSyn)
# sum(blnORFProper) / length(GRNonSyn)

# Check that most of first nucleotides match among non synonymous 
SynPos    <- seq(3, 21, 3)
NonSynPos <- setdiff(1:21, SynPos)
blnORF1_Nonsyn <- GRNonSynInt@elementMetadata@listData$mcols.ORFType == "ORF1"
blnORF2_Nonsyn <- GRNonSynInt@elementMetadata@listData$mcols.ORFType == "ORF2"
idxMismatchNonSyn <- lapply(NonSynPos, function(i){
  cat("\nAnalyzing sequence position", i, "\n")
  Nuc1 <- substr(StartSeqORF1_long, i, i)
  Nuc2 <- substr(StartSeqORF2_long, i, i)
  GR_Pos_ORF1 <- GRNonSynInt[blnORF1_Nonsyn &
                               GRNonSynInt@elementMetadata@listData$mcols.ORFPos0 == i - 1]
  Nuc_ORF1 <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, GR_Pos_ORF1))
  GR_Pos_ORF2 <- GRNonSynInt[blnORF2_Nonsyn &
                               GRNonSynInt@elementMetadata@listData$mcols.ORFPos0 == i - 1]
  Nuc_ORF2 <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, GR_Pos_ORF2))
  idxMismatch_ORF1 <- which(Nuc_ORF1 != Nuc1)
  cat(length(idxMismatch_ORF1), "out of", length(Nuc_ORF1), "L1s with mismatching nuc on ORF1\n")
  idxMismatch_ORF2 <- which(Nuc_ORF2 != Nuc2)
  cat(length(idxMismatch_ORF2), "out of", length(Nuc_ORF2), "L1s with mismatching nuc on ORF2\n")
  list(idxMismatch_ORF1 = idxMismatch_ORF1,
       idxMismatch_ORF2 = idxMismatch_ORF2)
  
})

blnORF1_Syn <- GRSynInt@elementMetadata@listData$mcols.ORFType == "ORF1"
blnORF2_Syn <- GRSynInt@elementMetadata@listData$mcols.ORFType == "ORF2"
idxMismatchSyn <- lapply(SynPos, function(i){
  cat("\nAnalyzing sequence position", i, "\n")
  Nuc1 <- substr(StartSeqORF1_long, i, i)
  Nuc2 <- substr(StartSeqORF2_long, i, i)
  GR_Pos_ORF1 <- GRSynInt[blnORF1_Syn &
                               GRSynInt@elementMetadata@listData$mcols.ORFPos0 == i - 1]
  Nuc_ORF1 <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, GR_Pos_ORF1))
  GR_Pos_ORF2 <- GRSynInt[blnORF2_Syn &
                               GRSynInt@elementMetadata@listData$mcols.ORFPos0 == i - 1]
  Nuc_ORF2 <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, GR_Pos_ORF2))
  idxMismatch_ORF1 <- which(Nuc_ORF1 != Nuc1)
  cat(length(idxMismatch_ORF1), "out of", length(Nuc_ORF1), "L1s with mismatching nuc on ORF1\n")
  idxMismatch_ORF2 <- which(Nuc_ORF2 != Nuc2)
  cat(length(idxMismatch_ORF2), "out of", length(Nuc_ORF1), "L1s with mismatching nuc on ORF2\n")
  list(idxMismatch_ORF1 = idxMismatch_ORF1,
       idxMismatch_ORF2 = idxMismatch_ORF2)
  
})

# Plot how often a particular position within L1 is counted among non-
# synymous positions
ORFPosCount <- table(GRNonSynInt@elementMetadata@listData$mcols.ORFPos)
plot(as.numeric(names(ORFPosCount)), ORFPosCount)
plot(table(GRNonSynInt@elementMetadata@listData$mcols.ORFPos[
  GRNonSynInt@elementMetadata@listData$mcols.ORFType == "ORF1"]
))
plot(table(GRNonSynInt@elementMetadata@listData$mcols.ORFPos[
  GRNonSynInt@elementMetadata@listData$mcols.ORFType == "ORF2"]
))
cat("done!\n")


# Get ranges of genes, exons' coding sequences, etc.
Exons <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
Genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
Cds <- cds(TxDb.Hsapiens.UCSC.hg19.knownGene)
Promoters <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream = 5*10^3)
CdsTx <- cdsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "tx")
idxOLCdsL1 <- which(overlapsAny(CdsTx, L1GR))
NonSynStart   <- NULL
SynStart      <- NULL
ChrVCode_Gene <- NULL

# Calculate the number of basepairs that differ between indels and the reference
L1Variants_Indel$bpDiff <- sapply(1:nrow(L1Variants_Indel), function(x) {
  NC1 <- nchar(L1Variants_Indel$REF[x])
  Split2 <- strsplit(L1Variants_Indel$ALT[x], ",")[[1]]
  NC1 - max(nchar(Split2))
})

# Indicator variable of whether indel is divisible by 3 
L1Variants_Indel$bln3 <- (L1Variants_Indel$bpDiff %% 3) == 0

# Indicator variable of whether indel overlaps with an ORF
#L1Variants_Indel$blnORF <- overlapsAny(L1VarIndelGR, ORF12_GR)

# GRanges for full-length L1s, their ORFs and UTRs
L1GRFull     <- L1GR[width(L1GR) >= 6000]

# Indicator variable for whether Indel is in full-length L1
L1Variants_Indel$blnFull <- overlapsAny(L1VarIndelGR, L1GRFull)
L1Variants_IndelFull  <- L1Variants_Indel[L1Variants_Indel$blnFull, ]
L1VarIndelGRFull      <- L1VarIndelGR[L1Variants_Indel$blnFull]


# Subset GRanges for high quality SNPs
L1VarGR_HighQual <- L1VarGR[L1Variants$QUAL >= 100]
table(L1Variants$FILTER)

# Count the number of variants per L1
L1VarCount <- countOverlaps(L1GR, L1VarGR)
L1VarCount_Left  <- countOverlaps(L1GR_left, L1VarGR_Left)
L1VarCount_Right <- countOverlaps(L1GR_right, L1VarGR_Right)
L1VarCount_Flank <- L1VarCount_Left + L1VarCount_Right
cor.test(L1VarCount_Left, L1VarCount_Right) 

# Get GRanges of SNPs that do not deviate from HWE and are 
GR_HWE <- L1VarGR[!overlapsAny(L1VarGR, HWEGR)]

# Overlaps between individual bp and L1
OL_bpL1_MEDel <- findOverlaps(L1Cover_GR, MEDel_GR)
GRNonSynInt$mcols.ORFType == "ORF1"
OL_bpL1      <- findOverlaps(L1Cover_GR, L1GR)
OL_bpL1Var   <- findOverlaps(L1Cover_GR, L1VarGR)
OL_bpHWE     <- findOverlaps(L1Cover_GR, GR_HWE)
GRNonSynInt1 <- GRNonSynInt[GRNonSynInt$mcols.ORFType == "ORF2"]
GRSynInt1    <- GRSynInt[GRSynInt$mcols.ORFType == "ORF2"]

# Calculate number of synonymous and non-synonymouys SNPs per L1
OL_FullL1nonSyn  <- findOverlaps(L1GR[width(L1GR) >= 6000], GRNonSynInt1)
width(GRNonSynInt1)
OL_FullL1Syn  <- findOverlaps(L1GR[width(L1GR) >= 6000],  GRSynInt1)
VarNonVarNonSyn <- sapply(unique(OL_FullL1nonSyn@from), function(x){
  idxGR  <- OL_FullL1nonSyn@to[OL_FullL1nonSyn@from == x]
  blnVar <- overlapsAny(GRNonSynInt1[idxGR], L1VarGR)
  c(NrVar = sum(blnVar), NonVar = sum(!blnVar))
})
VarNonVarSyn <- sapply(unique(OL_FullL1Syn@from), function(x){
  idxGR  <- OL_FullL1Syn@to[OL_FullL1Syn@from == x]
  blnVar <- overlapsAny(GRSynInt1[idxGR], L1VarGR)
  c(NrVar = sum(blnVar), NonVar = sum(!blnVar))
})
all(unique(OL_FullL1nonSyn@from) == unique(OL_FullL1Syn@from))
plot(jitter(VarNonVarSyn[1,]), jitter(VarNonVarNonSyn[1,]),
     xlab = "Number of synonymous SNPs", 
     ylab = "Number of non-synonymous SNPs")
lines(c(0, 100), c(0, 200))
x <-1 
idx10 <- which((VarNonVarSyn[1,] + VarNonVarNonSyn[1,]) >= 10)
PVals <- sapply(idx10, function(x){
  fisher.test(cbind(VarNonVarNonSyn[,x], VarNonVarSyn[,x]))$p.value
})
hist(PVals, breaks = seq(0, 1, 0.05))
p.adjust(PVals)

# Add columns with different info
L1CoverTable$blnSNP       <- overlapsAny(L1Cover_GR, L1VarGR)
L1CoverTable$blnSNPPacBio <- overlapsAny(L1Cover_GR, L1VarGRPacBio)
L1CoverTable$blnSNPHG002  <- overlapsAny(L1Cover_GR, L1VarGRHG002)
L1CoverTable$blnSNP_HighQual <- overlapsAny(L1Cover_GR, L1VarGR_HighQual)
L1CoverTable$blnSNPHWE <- overlapsAny(L1Cover_GR, GR_HWE)
L1CoverTable$blnSNP_both <- L1CoverTable$blnSNPHWE & L1CoverTable$blnSNP_HighQual
L1CoverTable$Dist2Edge <- NA
L1CoverTable$Dist2Edge[OL_bpL1@from] <- pmin(abs(start(L1Cover_GR)[OL_bpL1@from] - 
                                     start(L1GR)[OL_bpL1@to]),
                               abs(start(L1Cover_GR)[OL_bpL1@from] - 
                                     end(L1GR)[OL_bpL1@to]))
L1CoverTable$Exons     <- overlapsAny(L1Cover_GR, Exons)
L1CoverTable$Genes     <- overlapsAny(L1Cover_GR, Genes)
L1CoverTable$Promoters <- overlapsAny(L1Cover_GR, Promoters)
L1CoverTable$TFB       <- overlapsAny(L1Cover_GR, unlist(GRangesList(TFBGR)))
L1CoverTable$L1VarCount_Flank <- NA
L1CoverTable$L1VarCount_Flank[OL_bpL1@from] <- L1VarCount_Flank[OL_bpL1@to]
L1CoverTable$PropMismatch     <- NA
L1CoverTable$PropMismatch[OL_bpL1@from] <- PropMismatch[OL_bpL1@to]
L1CoverTable$L1Width          <- NA
L1CoverTable$L1Width[OL_bpL1@from] <- L1Width[OL_bpL1@to]
L1CoverTable$blnFull          <- NA
L1CoverTable$blnFull[OL_bpL1@from] <- blnFull[OL_bpL1@to]
L1CoverTable$AlleleFreq          <- NA
L1CoverTable$AlleleFreq[OL_bpL1Var@from] <- L1Variants$AlleleFreq[OL_bpL1Var@to]
L1CoverTable$Freq <- 2*length(GTCols)
L1CoverTable$Freq[OL_bpL1_MEDel@from] <- MEDelCall$Freq[OL_bpL1_MEDel@to]
L1CoverTable$Coding <- overlapsAny(L1Cover_GR, c(GRSynInt, GRNonSynInt))
L1CoverTable$NonSyn <- overlapsAny(L1Cover_GR, GRNonSynInt) 
L1CoverTable$Syn    <- overlapsAny(L1Cover_GR, GRSynInt) 
L1CoverTable$NonSyn_ORF1 <- overlapsAny(L1Cover_GR, 
                                       GRNonSynInt[GRNonSynInt@elementMetadata@listData$mcols.ORFType == "ORF1"]) 
L1CoverTable$Syn_ORF1      <- overlapsAny(L1Cover_GR, 
                                        GRSynInt[GRSynInt@elementMetadata@listData$mcols.ORFType == "ORF1"]) 
L1CoverTable$NonSyn_ORF2 <- overlapsAny(L1Cover_GR, 
                                       GRNonSynInt[GRNonSynInt@elementMetadata@listData$mcols.ORFType == "ORF2"]) 
L1CoverTable$Syn_ORF2     <- overlapsAny(L1Cover_GR, 
                                        GRSynInt[GRSynInt@elementMetadata@listData$mcols.ORFType == "ORF2"]) 

L1CoverTable$NonSyn_Proper <- overlapsAny(L1Cover_GR, GRNonSynInt[c(blnORFProperInt,
                                                                    blnORFProperInt)])
L1CoverTable$Syn_Proper    <- overlapsAny(L1Cover_GR, GRSynInt[blnORFProperInt])

######################################################
#                                                    #
#           Analyze LD                               #
#                                                    #
######################################################

# Read in files with linkage disequilibrium
LD_Chr1 <- read.delim("file:///D:/OneDrive - American University of Beirut/L1polymORF/Data/LD_chr1.hap.ld")
LD_Chr2 <- read.delim("file:///D:/OneDrive - American University of Beirut/L1polymORF/Data/LD_chr2.hap.ld")
LD_Chr3 <- read.delim("file:///D:/OneDrive - American University of Beirut/L1polymORF/Data/LD_chr3.hap.ld")
#LD_Chr <- LD_Chr2
LD_Chr <- rbind(LD_Chr1, LD_Chr2, LD_Chr3)
LD_Chr$chromosome <- paste0("chr", LD_Chr$CHR)
hist(LD_Chr$Dprime)
hist(LD_Chr$D, breaks = seq(-0.3, 0.3, 0.001), ylim = c(0, 200))

# Create genomic ranges for pairs of SNPs
LD_Chr_GR1 <- makeGRangesFromDataFrame(LD_Chr,
                                        seqnames.field = "chromosome",
                                        start.field = "POS1",
                                        end.field = "POS1")
LD_Chr_GR2 <- makeGRangesFromDataFrame(LD_Chr,
                                        seqnames.field = "chromosome",
                                        start.field = "POS2",
                                        end.field = "POS2")
# Add boolean variables indicating whether SNP is synonymous or 
# non-synynymous
LD_Chr$NonSyn1 <- overlapsAny(LD_Chr_GR1, GRNonSynInt) 
LD_Chr$NonSyn2 <- overlapsAny(LD_Chr_GR2, GRNonSynInt) 
LD_Chr$blnBothNonSyn <- LD_Chr$NonSyn1 & LD_Chr$NonSyn2

# Match SNPs to individual LINE-1s
OL_LD1 <- findOverlaps(LD_Chr_GR1, L1GR)
OL_LD2 <- findOverlaps(LD_Chr_GR2, L1GR)
LD_Chr$L1ID1 <- NA
LD_Chr$L1ID1[OL_LD1@from] <- OL_LD1@to
LD_Chr$L1ID2 <- NA
LD_Chr$L1ID2[OL_LD2@from] <- OL_LD2@to
blnSameL1 <- LD_Chr$L1ID1 == LD_Chr$L1ID2
sum(!blnSameL1, na.rm = T)
t.test(D ~ blnBothNonSyn, data = LD_Chr[blnSameL1, ])
t.test(D ~ blnBothNonSyn, data = LD_Chr[!blnSameL1, ])
wilcox.test(Dprime ~ blnBothNonSyn, data = LD_Chr[blnSameL1, ])
wilcox.test(Dprime ~ blnBothNonSyn, data = LD_Chr[!blnSameL1, ], FUN = mean)
aggregate(Dprime ~ blnBothNonSyn, data = LD_Chr[blnSameL1, ], FUN = mean)
aggregate(Dprime ~ blnBothNonSyn, data = LD_Chr[!blnSameL1, ], FUN = mean)
wilcox.test(D ~ blnBothNonSyn, data = LD_Chr[blnSameL1, ])
wilcox.test(D ~ blnBothNonSyn, data = LD_Chr[!blnSameL1, ], FUN = mean)
aggregate(D ~ blnBothNonSyn, data = LD_Chr[blnSameL1, ], FUN = mean)
aggregate(D ~ blnBothNonSyn, data = LD_Chr[!blnSameL1, ], FUN = mean)

hist(LD_Chr$Dprime)
hist(LD_Chr$D[LD_Chr$blnBothNonSyn], 
     breaks = seq(-0.3, 0.3, 0.001), ylim = c(0, 100))
hist(LD_Chr$D[!LD_Chr$blnBothNonSyn], 
     breaks = seq(-0.3, 0.3, 0.001), ylim = c(0, 2000))


######################################################
#                                                    #
#           Perform logistic regression              #
#           to compare 1000G and PacBio              #
#                                                    #
######################################################

# Get number of different types of SNPs
sum(L1CoverTable$blnSNP)
sum(L1CoverTable$blnSNPPacBio)
sum(L1CoverTable$blnSNP_HighQual)
sum(L1CoverTable$blnSNPHWE)
sum(L1CoverTable$blnSNP_both)
sum(L1CoverTable$blnSNPHG002)

# Check whether SNPs from different datasets are associated
SNPCrossTab <- table(L1CoverTable$blnSNPPacBio, L1CoverTable$blnSNP)
SNPCrossTab[2,]/colSums(SNPCrossTab)
SNPCrossTab_both <- table(L1CoverTable$blnSNPPacBio, L1CoverTable$blnSNP_both)
SNPCrossTab_both[2,]/colSums(SNPCrossTab_both)
SNPCrossTab_HG002 <- table(L1CoverTable$blnSNPPacBio, L1CoverTable$blnSNPHG002)
SNPCrossTab_HG002[2,]/colSums(SNPCrossTab_HG002)

# analyze whether the probability to have a SNP that is not PACBio depends
# on the interesting covariates
SNPLogReg_PacBio1 <- bigglm(blnSNPPacBio ~   TriNuc + L1VarCount_Flank + CoverMean +
                                         PropMismatch + Genes + Exons + Promoters + 
                                         blnFull + Syn + Coding + 
                                         Coding*blnFull + Syn*blnFull,
                                       data = L1CoverTable[L1CoverTable$blnSNPHG002 |
                                                             L1CoverTable$blnSNPPacBio, ],
                                       family = binomial(), chunksize = 3*10^4,
                                       maxit = 20)
summary(SNPLogReg_PacBio1)
SNPLogReg_PacBio1_glm <- glm(blnSNPPacBio ~   TriNuc + L1VarCount_Flank + CoverMean +
                              PropMismatch + Genes + Exons + Promoters + 
                              blnFull + NonSyn + Coding + 
                              Coding*blnFull + NonSyn*blnFull,
                            data = L1CoverTable[L1CoverTable$blnSNP_both |
                                                  L1CoverTable$blnSNPPacBio, ],
                            family = binomial())

summary(SNPLogReg_PacBio1_glm)
SNPLogReg_PacBio2 <- bigglm(blnSNPPacBio ~   CoverMean +
                                 NonSyn + Coding ,
                               data = L1CoverTable[L1CoverTable$blnSNP_both, ],
                               family = binomial(), chunksize = 3*10^4,
                               maxit = 20)
summary(SNPLogReg_PacBio2)
SNPLogReg_NonPacBio3 <- bigglm(blnSNPPacBio ~   TriNuc + L1VarCount_Flank + CoverMean +
                                 PropMismatch + Genes + Exons + Promoters + 
                                 blnFull + NonSyn + Coding + 
                                 Coding*blnFull + NonSyn*blnFull,
                               data = L1CoverTable[L1CoverTable$blnSNPHG002, ],
                               family = binomial(), chunksize = 3*10^4,
                               maxit = 20)
summary(SNPLogReg_NonPacBio3)

######################################################
#                                                    #
#           Perform logistic regression              #
#                                                    #
######################################################

# Perform analysis with interaction
cat("Performing regression analysis with all SNPs... ")
SNPLogRegInt <- bigglm(blnSNP_both ~  TriNuc + L1VarCount_Flank + CoverMean +
                         PropMismatch + Genes + Exons + Promoters + Dist2Edge +
                         blnFull + Syn + Coding + 
                         Coding*blnFull + Syn*blnFull,
                       data = L1CoverTable, 
                       family = binomial(), chunksize = 3*10^4,
                       maxit = 20)
summary(SNPLogRegInt)
cat("done!\n")

# Perform analysis with interaction and different coefficients per ORF
cat("Performing regression analysis with all SNPs and different coefficients per ORF... ")
SNPLogRegInt_byORF <- bigglm(blnSNP_both ~  TriNuc + L1VarCount_Flank + CoverMean +
                         PropMismatch + Genes + Exons + Promoters + #TFB +
                         blnFull + Syn_ORF1 + Syn_ORF2 + Coding + Freq +
                         Coding*blnFull + Syn_ORF1*blnFull + Syn_ORF2*blnFull,
                       data = L1CoverTable, 
                       family = binomial(), chunksize = 3*10^4,
                       maxit = 20)
summary(SNPLogRegInt_byORF)
cat("done!\n")

# Perform analysis with interaction with SNPs from PacBio genome
cat("Performing regression analysis with SNPs from PacBio genome... ")
SNPLogRegPacBio <- bigglm(blnSNPPacBio ~  TriNuc + L1VarCount_Flank + CoverMean +
                            PropMismatch + Genes + Exons + Promoters + Dist2Edge + 
                            blnFull + Syn + Coding + 
                            Coding*blnFull + Syn*blnFull,
                       data = L1CoverTable, 
                       family = binomial(), chunksize = 3*10^4,
                       maxit = 20)
summary(SNPLogRegPacBio)
cat("done!\n")

# cat("Performing regression analysis with coding sequences only ... \n")
SNPLogReg_CodeOnly <- bigglm(blnSNP_both ~  TriNuc + L1VarCount_Flank + CoverMean +
                               PropMismatch + Genes + Promoters +
                               blnFull + Syn + blnFull*Syn,
                             data = L1CoverTable[(L1CoverTable$Syn | L1CoverTable$NonSyn), ],
                             family = binomial(), chunksize = 3*10^4,
                             maxit = 20)
summary(SNPLogReg_CodeOnly)
SNPLogReg_CodeOnlyPacBio <- bigglm(blnSNPPacBio ~  TriNuc + L1VarCount_Flank + CoverMean +
                               PropMismatch + Genes + Promoters +
                                 blnFull + Syn + blnFull*Syn,
                               data = L1CoverTable[(L1CoverTable$Syn | L1CoverTable$NonSyn), ],
                             family = binomial(), chunksize = 3*10^4,
                             maxit = 20)
summary(SNPLogReg_CodeOnlyPacBio)


# cat("Performing regression analysis with coding sequences on full-length L1 only ... \n")
SNPLogReg_CodeOnlyFull <- bigglm(blnSNP_both ~  TriNuc + L1VarCount_Flank + CoverMean +
                         PropMismatch + Genes + Promoters + Syn,
                    data = L1CoverTable[L1CoverTable$blnFull & L1CoverTable$Coding, ],
                    family = binomial(), chunksize = 3*10^4,
                    maxit = 20)
summary(SNPLogReg_CodeOnlyFull)
SNPLogReg_CodeOnlyFullPacBio <- bigglm(blnSNPPacBio ~  TriNuc + L1VarCount_Flank + CoverMean +
                                 PropMismatch + Genes + Promoters +
                                   Syn,
                                data = L1CoverTable[L1CoverTable$blnFull & (L1CoverTable$Syn | L1CoverTable$NonSyn), ],
                                family = binomial(), chunksize = 3*10^4,
                                maxit = 20)
summary(SNPLogReg_CodeOnlyFullPacBio)

# cat("Performing regression analysis with coding sequences on fragement L1 only ... \n")
SNPLogReg_CodeOnlyNotFull <- bigglm(blnSNP_both ~  TriNuc + L1VarCount_Flank + CoverMean +
                               PropMismatch + Genes + Promoters +
                                 Syn,
                             data = L1CoverTable[!L1CoverTable$blnFull & (L1CoverTable$Syn | L1CoverTable$NonSyn), ],
                             family = binomial(), chunksize = 3*10^4,
                             maxit = 20)
summary(SNPLogReg_CodeOnlyNotFull)
SNPLogReg_CodeOnlyNotFullPacBio <- bigglm(blnSNPPacBio ~  TriNuc + L1VarCount_Flank + CoverMean +
                                     PropMismatch + Genes + Promoters +
                                       Syn,
                                   data = L1CoverTable[!L1CoverTable$blnFull & (L1CoverTable$Syn | L1CoverTable$NonSyn), ],
                                   family = binomial(), chunksize = 3*10^4,
                                   maxit = 20)
summary(SNPLogReg_CodeOnlyNotFullPacBio)

cat("Performing regression analysis with full L1 only ... \n")
SNPLogReg_Full <- bigglm(blnSNP_both ~  TriNuc + L1VarCount_Flank + CoverMean +
                           PropMismatch + Genes + Exons + Promoters + 
                              Coding + Syn, 
                                 data = L1CoverTable[L1CoverTable$blnFull, ], 
                                family = binomial(), chunksize = 3*10^4,
                                maxit = 20)
summary(SNPLogReg_Full)
SNPLogReg_FullbyORF <- bigglm(blnSNP_both ~  TriNuc + L1VarCount_Flank + CoverMean +
                                PropMismatch + Genes + Exons + Promoters + 
                              Coding + Syn_ORF1 + Syn_ORF2, 
                            data = L1CoverTable[L1CoverTable$blnFull, ], 
                            family = binomial(), chunksize = 3*10^4,
                            maxit = 20)
summary(SNPLogReg_FullbyORF)
SNPLogReg_FullPacBio <- bigglm(blnSNPPacBio ~  TriNuc + L1VarCount_Flank + CoverMean +
                                 PropMismatch + Genes + Exons + Promoters + 
                                 Coding + Syn, 
                                   data = L1CoverTable[L1CoverTable$blnFull, ], 
                                   family = binomial(), chunksize = 3*10^4,
                                   maxit = 20)
summary(SNPLogReg_FullPacBio)

cat("done!\n")

# Function to create summary dataframe
NZeroDigits <- function(x)  min(c(10, which(abs(x) %/% 10^-(1:10) > 0)))
SummaryDF <- function(LM, MinP = 10^-5){
  Summary <- summary(LM)
  SummaryDF <- as.data.frame(Summary$mat)
  NZeroDigitM <- sapply(1:nrow(SummaryDF), function(x){
    c(Coef = NZeroDigits(SummaryDF$Coef[x]),
      ExpCoef = NZeroDigits(exp(SummaryDF$Coef[x])),
      p = NZeroDigits(SummaryDF$p[x]))
  })
  SummaryDF$Predictor <- row.names(SummaryDF)
  SummaryDF$Coef      <- round(SummaryDF$Coef, NZeroDigitM["Coef",] + 1) 
  SummaryDF$ExpCoef   <- round(exp(SummaryDF$Coef), NZeroDigitM["ExpCoef",] + 1) 
  SummaryDF$p         <- round(SummaryDF$p, NZeroDigitM["p",] + 1) 
  SummaryDF$p         <- round(SummaryDF$p, NZeroDigitM["p",] + 1) 
  RepChar <- paste(c("<0.", rep(0, NZeroDigits(MinP) - 1), "1"), collapse = "")
  SummaryDF$p[SummaryDF$p < MinP] <- RepChar
  #  SummaryDF[,c("Predictor", "Coef", "ExpCoef", "p")]
  SummaryDF[,c("Predictor", "Coef", "p")]
}


# Export data frame with regression results
SNPLogRegInt_Summary   <- SummaryDF(SNPLogRegInt)
SNPLogReg_Full_Summary <- SummaryDF(SNPLogReg_CodeOnlyFull)
SNPLogReg_NotFull_Summary  <- SummaryDF(SNPLogReg_CodeOnlyNotFull)
SNPLogRegIntPacBio_Summary   <- SummaryDF(SNPLogRegPacBio)
SNPLogReg_FullPacBio_Summary <- SummaryDF(SNPLogReg_CodeOnlyFullPacBio)
SNPLogReg_NotFullPacBio_Summary  <- SummaryDF(SNPLogReg_CodeOnlyNotFullPacBio)
SNPLogRegMergedAll <- merge(SNPLogRegInt_Summary, 
                         SNPLogRegIntPacBio_Summary,
                         by = "Predictor")
SNPLogRegMergedFull <- merge(SNPLogReg_Full_Summary, 
                             SNPLogReg_FullPacBio_Summary,
                         by = "Predictor", all = T)
SNPLogRegMergedFragm <- merge(SNPLogReg_NotFull_Summary, 
                              SNPLogReg_NotFullPacBio_Summary,
                         by = "Predictor", all = T)
SNPLogRegMerged <- merge(SNPLogRegMergedAll, 
                         SNPLogRegMergedFull,
                         by = "Predictor", all = T)
SNPLogRegMerged <- merge(SNPLogRegMerged, 
                         SNPLogRegMergedFragm,
                         by = "Predictor", all = T)

cat("Writing regression results to", ResultPathCombined, "\n")
write.csv(SNPLogRegMerged, ResultPathCombined)
write.csv(SNPLogRegMergedAll, ResultPathAll)
write.csv(SNPLogRegMergedFull, ResultPathFull)


######################################################
#                                                    #
#     Plot SNP probability along full-length L1      #
#                                                    #
######################################################


# Analyze the proportion of SNPs on non-synonymous sites 
idxFull <- which(width(L1GR) >= 6000)
NonsynEffect <- sapply(idxFull, function(i){
  L1Subset <- L1CoverTable[OL_bpL1@from[OL_bpL1@to == i], ]
  L1Subset <- L1Subset[L1Subset$Coding, ]
  if(any(L1Subset$blnSNP)){
    print(i)
    LReg <- glm(blnSNP ~ NonSyn, data = L1Subset, family = binomial)
    LRegCoeffs <- summary(LReg)$coefficients
    LRegCoeffs <- rbind(LRegCoeffs, rep(NA, 4))
    LRegCoeffs[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
  } else {
    c(NA, NA, NA)
  }
})
hist(NonsynEffect[1,])
hist(NonsynEffect[3,], breaks = seq(0, 1, 0.01))
PAdj <- p.adjust(NonsynEffect[3,])
sum(PAdj < 0.05, na.rm = T)
min(PAdj, na.rm = T)
NonsynEffect[,which.min(PAdj)]

# Get mean number of SNPs per ORF type position
AggPerORFPos <- AggDataFrame(DF = L1CoverTable[L1CoverTable$NonSyn | L1CoverTable$Syn,], 
                            GroupCol = c("blnFull", "NonSyn"), 
                            MeanCols = c("blnSNP", "blnSNP_both", "blnSNPPacBio",
                                         "blnSNPHG002"), 
                            LengthCols = "blnSNP", 
                            VarCols = c("blnSNP", "blnSNP_both", "blnSNPPacBio",
                                        "blnSNPHG002"))

# Get mean number of SNPs per full-length and fragment L1
AggPerL1Pos_FullFrag <- AggDataFrame(DF = L1CoverTable, 
                             GroupCol = "blnFull", 
                             MeanCols = c("blnSNP", "blnSNP_both", "blnSNPPacBio",
                                          "blnSNPHG002"), 
                             LengthCols = "blnSNP", 
                             VarCols = c("blnSNP", "blnSNP_both", "blnSNPPacBio",
                                         "blnSNPHG002"))

# Get mean number of SNPs per full-length and fragment L1
AggPerL1Pos_ORFvsUTR <- AggDataFrame(DF = L1CoverTable[L1CoverTable$blnFull, ], 
                                GroupCol = "Coding", 
                                MeanCols = c("blnSNP", "blnSNP_both", "blnSNPPacBio",
                                             "blnSNPHG002"), 
                                LengthCols = "blnSNP", 
                                VarCols = c("blnSNP", "blnSNP_both", "blnSNPPacBio",
                                            "blnSNPHG002"))

# Function to create barplot of mean SNPs and error bars
PlotMeanSNP <- function(AggDF, NameV, Col2Plot  = "blnSNP",
                        YLab = "", YLim = c(0, 0.02), Main = "",
                        PlotP = NULL, Border = T){
  MeanCol <- paste0(Col2Plot, "_mean")
  VarCol  <- paste0(Col2Plot, "_var")
  BP <- barplot(AggDF[,MeanCol], main = Main,
                ylab = YLab,
                names = NameV, ylim = YLim, border = Border) 
  AddErrorBars(Vertical = T, MidX = BP, 
               MidY = AggDF[,MeanCol],
               ErrorRange = sqrt(AggDF[,VarCol]/ AggDF$blnSNP_N),
               TipWidth = 0.1)
  if (!is.null(PlotP)){
    LX <- 1.1 * max(AggDF[,MeanCol], na.rm = T)
    TX <- 1.2 * max(AggDF[,MeanCol], na.rm = T)
    lines(BP, c(LX, LX))
    text(x = mean(BP), y = TX, PlotP, cex = 1)
  }
  
  return(BP)
}

######################################################
#                                                    #
#     Plot SNP probability along full-length L1      #
#     With SNPs and coverage
#                                                    #
######################################################

# Get sequences and create a character of sequence names
L1Seq    <- getSeq(BSgenome.Hsapiens.UCSC.hg19, L1GR)
SeqNames <- paste(as.vector(seqnames(L1GR)), start(L1GR), end(L1GR), sep = "_")

# Form different subsets and write them out as fasta files
L1Aligned <- read.fasta(file = "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1seqHg19_minLength6000_aligned_withConsens.fas")

# Create a distance matrix for alignment
L1AlDist  <- dist.alignment(read.alignment(file = "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1seqHg19_minLength6000_aligned_withConsens.fas",
                                           format = "fasta"))
L1AlDist <- as.matrix(L1AlDist)
diag(L1AlDist) <- NA

# Get indicies of sequenced positions
idxSeqPosList <- lapply(L1Aligned, function(x) which(x != "-"))

# Identify start and end of ORF1 and ORF2 on alignment
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
L1GRFull <- L1GR[L1Match[-length(L1Match)]]

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
                        idxInGR    = OL_Pos@to[blnL1],
                        idxL1 = x)
    idxPosDf <- rbind(idxPosDf, NewDf)
  }
  idxPosDf
}

# Get positions of SNPs in alignment
SNPposDF    <- GenPos2AlignPos(L1VarGR)
SNPposAlign <- SNPposDF$PosInAlign

# Check that nucleotides are the same
blnSameNuc <- sapply(1:1000, function(x){
  GenomeNuc <- getSeq(BSgenome.Hsapiens.UCSC.hg19, L1VarGR[SNPposDF$idxInGR[i]])[[1]]
  as.character(GenomeNuc)== 
    toupper(L1Aligned[[SNPposDF$idxL1[i]]][SNPposDF$PosInAlign[i]])
})
all(blnSameNuc)
                                   
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
PropIndel  <- colMeans(L1IndelMat)

FigDim = 4000
jpeg(filename = 'D:/L1ManuscriptFigures/NewFig1.jpg',
     width = FigDim, height = FigDim, pointsize = FigDim/480*12,
     quality = 100)
layout(matrix(c(1, 1, 2, 3, 4, 5), 3, 2, byrow = TRUE))
plot(c(-300, max(SNPposAlign)), c(0, 1), type= "n", xlab = "", ylab="",xaxt="n",frame=F,
     yaxt = "n", main = "a")
segments(1, 0.05, max(SNPposAlign), 0.05) # UTRs
rect(c(ORF1Start, ORF2Start), c(0, 0), c(ORF1End, ORF2End), 0.1, 
     border = "black", col ="lightgrey") # ORFs
text(0.5 * c(ORF1Start + ORF1End, ORF2Start + ORF2End), 0.05, c("ORF1", "ORF2"), cex = 0.75)

segments(SNPposAlign, 0.15, SNPposAlign, 0.25, col= rgb(0, 0, 0, 0.1), lwd = 0.4)
text(c(-300, -300), c(0.2, 0.35), c("SNP", "Indel"), cex = 1)
ncol(L1IndelMat)
segments(1:ncol(L1IndelMat), 0.3, 1:ncol(L1IndelMat), 0.4, col= rgb(0, 0, 0, 0.3*(1 - PropIndel)), 
         lwd = 0.1)
lines(CoverMeanPerAlign$Group.1[CoverMeanPerAlign$Count >= 10], 
      0.1 + CoverMeanPerAlign$x[CoverMeanPerAlign$Count >= 10] / 10)
min(CoverMeanPerAlign$x[CoverMeanPerAlign$Count >= 10])
max(CoverMeanPerAlign$x[CoverMeanPerAlign$Count >= 10])
axis(2, at = 0.1 + seq(0.4, 0.8, 0.2), labels = seq(4, 8, 2))
mtext("Coverage", side = 2, line = 2,  at = 0.7)

# Plot SNP density for different L1 regions
#layout(rbind(c(1, 2), c(3, 3)))
par(page = F, mai = c(0.5, 1, 0.2, 0.1) * FigDim/480)

# Get mean number of SNPs per full-length and fragment L1
PlotMeanSNP(AggPerL1Pos_FullFrag, Col2Plot  = "blnSNPPacBio",
            NameV = c("Fragment", "Full-length"),
            Main = "b", YLim = c(0, 0.003), YLab = "SNPs per LINE-1 bp", Border = NA,
            PlotP = "P < 0.0001")
PlotMeanSNP(AggPerL1Pos_ORFvsUTR, NameV = c("UTR", "ORF"),
            Main = "c", YLim = c(0, 0.025), Border = NA, PlotP = "P < 0.0001")
PlotMeanSNP(AggPerORFPos[!AggPerORFPos$blnFull, ], 
            Col2Plot  = "blnSNPPacBio",
            NameV = c("synonymous", "non-synonymous"),
            YLim = c(0, 0.003), Main = "d", YLab = "SNPs per LINE-1 bp", Border = NA)
PlotMeanSNP(AggPerORFPos[AggPerORFPos$blnFull, ], 
            Col2Plot  = "blnSNPPacBio",
            NameV = c("synonymous", "non-synonymous"),
            YLim = c(0, 0.003), Main = "e", Border = NA,
            PlotP = "P = 0.003")
dev.off()


######################################################
#                                                    #
#     Plot SNP probability along full-length L1      #
#     for 1000 genome and PacBio data
#                                                    #
######################################################

# Get positions of PacBio SNPs in alignment
SNPposDF_PacBio    <- GenPos2AlignPos(L1VarGRPacBio)
SNPposAlign_PacBio <- SNPposDF_PacBio$PosInAlign
SNPposDF_HG002     <- GenPos2AlignPos(L1VarGRHG002)
SNPposAlign_HG002  <- SNPposDF_HG002$PosInAlign

# Calculate the number of SNPs per position
SmoothedSNPFreq <- function(SNPpos){
  SNPCount <- rep(0, length(PropIndel))
  SNPTable <- table(SNPpos)
  SNPCount[as.numeric(names(SNPTable))] <- SNPTable
  RelSNPCount <- SNPCount / (1 - PropIndel) / length(L1Aligned)
  supsmu(x = 1:length(SNPCount), y = RelSNPCount)
}
RelSNPCountSmoothed <- SmoothedSNPFreq(SNPposAlign)
RelSNPCountSmoothed_PacBio <- SmoothedSNPFreq(SNPposAlign_PacBio)
RelSNPCountSmoothed_HG002 <- SmoothedSNPFreq(SNPposAlign_HG002)

# Set up figure
FigDim = 4000
jpeg(filename = 'D:/L1ManuscriptFigures/NewNewFig1.jpg',
     width = FigDim, height = FigDim, pointsize = FigDim/480*12,
     quality = 100)
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE))
par(mai = c(0.5, 1, 0.2, 0.1) * FigDim/480)
plot(c(-300, max(SNPposAlign)), c(0, 0.06), type= "n", xlab = "", ylab="",
     xaxt="n",frame = F,
     yaxt = "n", main = "a")

segments(1, 0.005, max(SNPposAlign), 0.005) # UTRs
rect(c(ORF1Start, ORF2Start), c(0, 0), c(ORF1End, ORF2End), 0.01, 
     border = "black", col ="lightgrey") # ORFs
text(0.5 * c(ORF1Start + ORF1End, ORF2Start + ORF2End), 0.005, c("ORF1", "ORF2"), cex = 0.75)

lines(RelSNPCountSmoothed$x, RelSNPCountSmoothed$y + 0.02,
      lwd = FigDim/480)
lines(RelSNPCountSmoothed_PacBio$x, RelSNPCountSmoothed_PacBio$y + 0.02,
      lty =2, lwd = FigDim/480)
axis(2, at = seq(0.02, 0.06, 0.02), labels = seq(0, 0.04, 0.02))
mtext("SNPs per LINE-1 bp", side = 2, line = 3,  at = 0.04)

# Plot SNP density for different L1 regions
# par(page = F, mai = c(0.5, 1, 0.2, 0.1) * FigDim/480)
#par(page = F, mai = c(0.5, 1, 0.2, 0.1) )

# Get mean number of SNPs per full-length and fragment L1
PlotMeanSNP(AggPerORFPos[AggPerORFPos$blnFull, ], 
            Col2Plot  = "blnSNPPacBio",
            NameV = c("synonymous", "non-synonymous"),
            YLim = c(0, 0.003), Main = "b", 
            YLab = "SNPs per LINE-1 bp", Border = NA,
            PlotP = "P = 0.001")
PlotMeanSNP(AggPerORFPos[AggPerORFPos$blnFull, ], 
            Col2Plot  = "blnSNP_both",
            NameV = c("synonymous", "non-synonymous"),
            YLim = c(0, 0.003), Main = "c", Border = NA,
            PlotP = "P = 0.00001")
dev.off()


save.image("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1VariantCount.RData")
load("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1VariantCount.RData")
# 

######################################################
#                                                    #
#           Comparing SNPs with alignment            #
#                                                    #
######################################################

# Susbet to obtain SNPs with measured allele frequencies
blnSubset <- L1CoverTable$blnFull & 
  (!is.na(L1CoverTable$AlleleFreq)) & (L1CoverTable$NonSyn |L1CoverTable$Syn)
L1Cover_GRsubset <- L1Cover_GR[blnSubset]

# Get variation of SNP positions in L1 alignment
AlPosSNPsubset <- GenPos2AlignPos(L1Cover_GRsubset)

# Nucleotide on consensus sequence
ConsensNuc <- toupper(L1Aligned$L1HS_L1_Homo_sapiens[AlPosSNPsubset$PosInAlign])
AlPosSNPsubset$PosInAlign
AlPosSNPsubset$idxL1
x <- 1
ClosestNuc <- sapply(1:nrow(AlPosSNPsubset), function(x){
  idxL1       <- AlPosSNPsubset$idxL1[x]
  CurrentDist <- L1AlDist[idxL1,]
  idxClostest <- which.min(CurrentDist)
  toupper(L1Aligned[[idxClostest]][AlPosSNPsubset$PosInAlign[x]])
})
# Get indices of nonsynonymous SNPs where the alternative nucleotide is equal
# to the nucleotide on the L1 consensus
OL_L1subsetPacBio   <- findOverlaps(L1Cover_GRsubset, L1VarGRPacBio)
idxAltConsensPacBio <- OL_L1subsetPacBio@from[ConsensNuc[OL_L1subsetPacBio@from] == 
                         L1VarPacBio$ALT[OL_L1subsetPacBio@to]] 
idxAltClosestPacBio <- OL_L1subsetPacBio@from[ClosestNuc[OL_L1subsetPacBio@from] == 
                                                L1VarPacBio$ALT[OL_L1subsetPacBio@to]] 
OL_L1subset   <- findOverlaps(L1Cover_GRsubset, L1VarGR)
idxAltConsens <- OL_L1subset@from[ConsensNuc[OL_L1subset@from] == 
                                  L1Variants$ALT[OL_L1subset@to]] 


######################################################
#                                                    #
#          Analyze allele frequencies                #
#                                                    #
######################################################

hist(L1CoverTable$AlleleFreq)
sum(L1CoverTable$AlleleFreq > 0.5, na.rm = T)

# Susbet to obtain SNPs with measured allele frequencies
L1CoverTableSubset      <- L1CoverTable[blnSubset, ]
L1CoverTableSubset$L1ID <- OL_bpL1@to[blnSubset]
L1CoverTableSubset$SNPSyn    <- L1CoverTableSubset$blnSNP & (!L1CoverTableSubset$NonSyn)
L1CoverTableSubset$SNPNonSyn <- L1CoverTableSubset$blnSNP & L1CoverTableSubset$NonSyn
L1CoverTableSubset$SNPSynPacBio    <- L1CoverTableSubset$blnSNPPacBio & (!L1CoverTableSubset$NonSyn)
L1CoverTableSubset$SNPNonSynPacBio <- L1CoverTableSubset$blnSNPPacBio & L1CoverTableSubset$NonSyn

# Replace SNP frequency by the frequency of the reference allele when the
# alternative allele is equal to consensus (i.e. likely to be ancestral)
L1CoverTableSubset$AFDerived <- L1CoverTableSubset$AlleleFreq
L1CoverTableSubset$AFDerived[idxAltClosestPacBio] <- 
  (1 - L1CoverTableSubset$AlleleFreq)[idxAltClosestPacBio]
# L1CoverTableSubset$AFDerived[idxAltConsens] <- 
#   (1 - L1CoverTableSubset$AlleleFreq)[idxAltConsens]
hist(L1CoverTableSubset$AFDerived)
hist(log(L1CoverTableSubset$AFDerived))

# Regression 
LMAleleFreq <- lm(AFDerived ~ PropMismatch + Syn, data = L1CoverTableSubset)
summary(LMAleleFreq)

# Calculate mean allele frequency per L1 and position type
AlleleFreqPerL1andPosType <- aggregate(L1CoverTable$AlleleFreq[blnSubset],
                                       by =list(OL_bpL1@to[blnSubset], L1CoverTable$NonSyn[blnSubset]), 
                                       FUN = mean)

# Calculate number of synonymous and non-synonymous SNPs per L1
NSNPPerL1andPosType <- AggDataFrame(L1CoverTableSubset, GroupCol = "L1ID", MeanCols = "AlleleFreq", 
                                    SumCols = c("SNPSyn", "SNPNonSyn", "SNPSynPacBio", "SNPNonSynPacBio"))
plot(NSNPPerL1andPosType$SNPSyn_sum, NSNPPerL1andPosType$SNPNonSyn_sum,
     col = rgb(0, 0, 0, alpha = 0.15), pch = 16,
     xlab = "Number of synonymous SNPs", ylab = "Number of non-synonymous SNPs")  
lines(c(0, 100), c(0, 200))
OL_L1SNO_L1Cat <- findOverlaps(L1GR[NSNPPerL1andPosType$L1ID],
                               L1CatalogGR_hg19)
NSNPPerL1andPosType$Activity <- NA
NSNPPerL1andPosType$Activity[OL_L1SNO_L1Cat@from] <- 
  L1CatalogL1Mapped$ActivityNum[OL_L1SNO_L1Cat@to]
idxHigh <- which(NSNPPerL1andPosType$Activity > 0)
points(NSNPPerL1andPosType$SNPSyn_sum[idxHigh], NSNPPerL1andPosType$SNPNonSyn_sum[idxHigh],
       col = "red")

# Plot L1 activity against non-synonymous SNPS
plot(NSNPPerL1andPosType$SNPNonSyn_sum, NSNPPerL1andPosType$Activity)
plot(NSNPPerL1andPosType$SNPNonSyn_sum + 
       NSNPPerL1andPosType$SNPSyn_sum, NSNPPerL1andPosType$Activity)

plot(jitter(NSNPPerL1andPosType$SNPSynPacBio_sum), 
     jitter(NSNPPerL1andPosType$SNPNonSynPacBio_sum),
     col = rgb(0, 0, 0, alpha = 0.15), pch = 16,
     xlab = "Number of synonymous SNPs", ylab = "Number of non-synonymous SNPs")  
lines(c(0, 100), c(0, 200))
L1GR[NSNPPerL1andPosType$L1ID[NSNPPerL1andPosType$SNPSyn_sum >= 15]]

load("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1CatalogGRanges.RData")
findOverlaps(L1GR[NSNPPerL1andPosType$L1ID[NSNPPerL1andPosType$SNPSyn_sum >= 20]],
             L1CatalogGR_hg19)
NSNPPerL1andPosType

# Observed difference in mean allele frequencies
L1CoverTableSubset <- L1CoverTableSubset[L1CoverTableSubset$blnSNPPacBio,]
ObsMeanDiff <- mean(L1CoverTableSubset$AFDerived[!L1CoverTableSubset$NonSyn]) -  
  mean(L1CoverTableSubset$AFDerived[L1CoverTableSubset$NonSyn]) 

# Sample mean differences
TotNonSyn <- sum(L1CoverTableSubset$NonSyn)
TotNuc    <- nrow(L1CoverTableSubset)
idxVect <- 1:TotNuc
NSamples <- 10000
SampledMeanDiffs <- sapply(1:NSamples, function(x){
  idxNonSyn <- sample.int(TotNuc, size = TotNonSyn)
  idxSyn    <- setdiff(idxVect, idxNonSyn)
  mean(L1CoverTableSubset$AFDerived[idxSyn]) -  
    mean(L1CoverTableSubset$AFDerived[idxNonSyn]) 
})
mean(SampledMeanDiffs >= ObsMeanDiff)
hist(SampledMeanDiffs)
# Observed difference in mean allele frequencies
ObsMedianDiff <- median(L1CoverTableSubset$AFDerived[!L1CoverTableSubset$NonSyn]) -  
  median(L1CoverTableSubset$AFDerived[L1CoverTableSubset$NonSyn]) 

# Sample mean differences
SampledMedianDiffs <- sapply(1:NSamples, function(x){
  idxNonSyn <- sample.int(TotNuc, size = TotNonSyn)
  idxSyn    <- setdiff(idxVect, idxNonSyn)
  median(L1CoverTableSubset$AFDerived[idxSyn]) -  
    median(L1CoverTableSubset$AFDerived[idxNonSyn]) 
})
sum(SampledMedianDiffs >= ObsMedianDiff) / NSamples

# Get 
L1Count <- table(AlleleFreqPerL1andPosType$Group.1)
L1Both  <- as.numeric(names(L1Count)[L1Count == 2])
SNPSPerL1 <- t(sapply(L1Both, function(x){
  blnL1 <- AlleleFreqPerL1andPosType$Group.1 == x
  c(NonSynFreq = AlleleFreqPerL1andPosType$x[AlleleFreqPerL1andPosType$Group.2 & blnL1],
    SynFreq = AlleleFreqPerL1andPosType$x[!AlleleFreqPerL1andPosType$Group.2 & blnL1])
}))

aggregate(L1CoverTable$blnSNP[blnSubset],
          by =list(OL_bpL1@to[blnSubset]), 
          FUN = sum)
NSNPPerL1 <- aggregate(L1CoverTable$blnSNP[OL_bpL1@from],
                       by =list(OL_bpL1@to), 
                       FUN = sum)
hist(NSNPPerL1$x, breaks = 0:500)
mean(NSNPPerL1$x)
sqrt(var(NSNPPerL1$x))


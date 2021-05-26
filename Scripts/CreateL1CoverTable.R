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
ResultPathCombined <- "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1VariantRegrResults_combined_2021-01-02.csv"
ResultPathAll      <- "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1VariantRegrResults_all_2021-01-02.csv"
ResultPathFull     <- "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1VariantRegrResults_full_2021-01-02.csv"
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

# Get the allele frequency
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

# Add a column that says which allele is ancestral
L1Variants$Ancestral <- NA
for (i in 1:nrow(L1Variants)){
  InfoSplit <- strsplit(L1Variants$INFO[i], ";")[[1]]
  AF <- grep("AF=", InfoSplit, value = T)
  AF <- AF[-grep("_AF=", AF)]
  AF <- as.numeric(strsplit(AF, "AF=")[[1]][2])
  AncAl  <- substr(grep("AA=", InfoSplit, value = T), 4, 4)
  if (length(AncAl) > 0){
    RefAlt <- L1Variants[i, c("REF", "ALT")]
    L1Variants$Ancestral[i] <- colnames(RefAlt)[match(AncAl, RefAlt)]
    
  }
}

# Add derived allele frequency whenever they are known
L1Variants$AlleleFreqDerived <- NA
idxRef <- which(L1Variants$Ancestral == "REF")
idxAlt <- which(L1Variants$Ancestral == "ALT")
L1Variants$AlleleFreqDerived[idxRef] <- L1Variants$AlleleFreq[idxRef]
L1Variants$AlleleFreqDerived[idxAlt] <- 1 - L1Variants$AlleleFreq[idxAlt]
hist(L1Variants$AlleleFreqDerived)

# Get allele frequency per population
AllFreqPerPop <- data.frame(t(sapply(L1Variants$INFO, function(x){
  InfoSplit <- strsplit(x, ";")[[1]]
  AFPop <- grep("_AF=", InfoSplit, value = T)
  AFPopNameVals <- sapply(AFPop, function(y) strsplit(y, "AF=")[[1]])
  AFPopVals <- as.numeric(AFPopNameVals[2,])
  names(AFPopVals) <- AFPopNameVals[1,]
  AFPopVals
})))


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
L1VarCount       <- countOverlaps(L1GR, L1VarGR)

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
length(AllFreqPerPop$EAS_[OL_bpL1Var@to])
length(L1CoverTable$blnSNP_EAS[OL_bpL1Var@from])
L1CoverTable$blnSNP_EAS  <- NA
L1CoverTable$blnSNP_AMR  <- NA
L1CoverTable$blnSNP_AFR  <- NA
L1CoverTable$blnSNP_EUR  <- NA
L1CoverTable$blnSNP_SAS  <- NA

L1CoverTable$blnSNP_EAS[OL_bpL1Var@from]  <- AllFreqPerPop$EAS_[OL_bpL1Var@to] > 0
L1CoverTable$blnSNP_AMR[OL_bpL1Var@from]  <- AllFreqPerPop$AMR_[OL_bpL1Var@to] > 0
L1CoverTable$blnSNP_AFR[OL_bpL1Var@from]  <- AllFreqPerPop$AFR_[OL_bpL1Var@to] > 0
L1CoverTable$blnSNP_EUR[OL_bpL1Var@from]  <- AllFreqPerPop$EUR_[OL_bpL1Var@to] > 0
L1CoverTable$blnSNP_SAS[OL_bpL1Var@from]  <- AllFreqPerPop$SAS_[OL_bpL1Var@to] > 0
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
L1CoverTable$AFDerived          <- NA
L1CoverTable$AFDerived[OL_bpL1Var@from] <- L1Variants$AlleleFreqDerived[OL_bpL1Var@to]
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
L1CoverTable$L1ID <- NA
L1CoverTable$L1ID[OL_bpL1@from] <- OL_bpL1@to
L1CoverTable$AlleleFreqSyn <- NA
L1CoverTable$AlleleFreqSyn[L1CoverTable$Syn] <- L1CoverTable$AlleleFreq[L1CoverTable$Syn] 
L1CoverTable$AlleleFreqNonSyn <- NA
L1CoverTable$AlleleFreqNonSyn[L1CoverTable$NonSyn] <- L1CoverTable$AlleleFreq[L1CoverTable$NonSyn] 

###############################################################
#                                                             #
#      Determine which non-synonymous SNPs change AA seq      #
#                                                             #
###############################################################

# Plot how often a particular position within L1 is counted among non-
# synymous positions
ORFPosCount1_NonSyn <- table(GRNonSynInt@elementMetadata@listData$mcols.ORFPos0[
  GRNonSynInt@elementMetadata@listData$mcols.ORFType == "ORF1"
])
ORFPosCount1_Syn <- table(GRSynInt@elementMetadata@listData$mcols.ORFPos0[
  GRSynInt@elementMetadata@listData$mcols.ORFType == "ORF1"
  ])

plot(as.numeric(names(ORFPosCount1_NonSyn)), ORFPosCount1_NonSyn, 
     xlab = "Position on L1", ylab = "Count", type = "l")
points(as.numeric(names(ORFPosCount1_Syn)), ORFPosCount1_Syn)

ORFPosCount2_NonSyn <- table(GRNonSynInt@elementMetadata@listData$mcols.ORFPos0[
  GRNonSynInt@elementMetadata@listData$mcols.ORFType == "ORF2"
  ])
ORFPosCount2_Syn <- table(GRSynInt@elementMetadata@listData$mcols.ORFPos0[
  GRSynInt@elementMetadata@listData$mcols.ORFType == "ORF2"
  ])

plot(as.numeric(names(ORFPosCount2_NonSyn)), ORFPosCount2_NonSyn, 
     xlab = "Position on L1", ylab = "Count", type = "l")
points(as.numeric(names(ORFPosCount2_Syn)), ORFPosCount2_Syn)

# Function to get genomic ranges of ORFs
GetORFGRs <- function(ORFType, MaxSize){
  
  ORFStartGR <- GRNonSynInt[GRNonSynInt@elementMetadata@listData$mcols.ORFPos0 == 0 &
                               GRNonSynInt@elementMetadata@listData$mcols.ORFType == ORFType]
  blnPosORF <- as.vector(strand(ORFStartGR)) == "+"
  ORFStartEndPos <- start(ORFStartGR) + MaxSize*blnPosORF - MaxSize*(!blnPosORF)
  ORF_GR <- GRanges(seqnames = seqnames(ORFStartGR), 
                     ranges = IRanges(start = pmin(start(ORFStartGR), ORFStartEndPos),
                                      end = pmax(start(ORFStartGR), ORFStartEndPos)),
                     strand = strand(ORFStartGR))
  ORFStrands <- as.vector(strand(ORF_GR))
  for(i in 1:length(ORF_GR)){
    GRSubset <- subsetByOverlaps(GRSynInt, ORF_GR[i])
    if (ORFStrands[i] == "+"){
      end(ORF_GR[i]) <- max(start(GRSubset))
      
    } else {
      start(ORF_GR[i]) <- min(start(GRSubset)) - 2
    }
  }
  ORF_GR
}

# Get genomic ranges
ORF1_GR <- GetORFGRs(ORFType = "ORF1", MaxSize = 1030)
ORF2_GR <- GetORFGRs(ORFType = "ORF2", MaxSize = 4000)

ORF1_Seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, ORF1_GR)
ORF2_Seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, ORF2_GR)
ORF1_AA  <- Biostrings::translate(ORF1_Seq)
ORF2_AA  <- Biostrings::translate(ORF2_Seq)

GRNonSyn_ORF1 <- GRNonSynInt[
  GRNonSynInt@elementMetadata@listData$mcols.ORFType == "ORF1"]

OL_NonSyn_ORF1   <- findOverlaps(GRNonSyn_ORF1, ORF1_GR)
OL_NonSyn_L1Var  <- findOverlaps(GRNonSyn_ORF1, L1VarGR)
OL_L1Var_ORF1    <- findOverlaps(L1VarGR, ORF1_GR)
idxORF1NonSynSNP <- which(L1CoverTable$blnSNP & L1CoverTable$NonSyn_ORF1)
OL_NonSyn_ORF1   <- findOverlaps(L1Cover_GR[idxORF1NonSynSNP], ORF1_GR)
StrandV_ORF1 <- as.vector(strand(ORF1_GR))[OL_NonSyn_ORF1@to]
PosOnORF1      <- start(L1Cover_GR[idxORF1NonSynSNP])[OL_NonSyn_ORF1@from] - 
                 start(ORF1_GR)[OL_NonSyn_ORF1@to] + 1
PosOnORF1[StrandV_ORF1 == "-"] <- (end(ORF1_GR)[OL_NonSyn_ORF1@to] -
  end(L1Cover_GR[idxORF1NonSynSNP])[OL_NonSyn_ORF1@from] + 1)[StrandV_ORF1 == "-"] 
min(PosOnORF1)
subseq(ORF1_Seq[OL_NonSyn_ORF1@to], start = PosOnORF1, end = PosOnORF1)
idxMatch <- match(idxORF1NonSynSNP, OL_bpL1Var@from)
L1Variants$ALT[OL_bpL1Var@to[idxMatch]]
L1Variants$REF[OL_bpL1Var@to[idxMatch[OL_NonSyn_ORF1@from]]]

startL1Cover_GR[idxORF1NonSynSNP]
ALTNuc <- L1Variants$ALT[]
OL_bpL1Var
L1Variants$REF
L1Variants$ALT
# Save image
#save.image("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1CoverTable.RData")
# The script below creates various summaries per genomic range

##########################################
#                                        #
#     Load packages                      #
#                                        #
##########################################

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)

##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

# Specify file paths
DataPath        <- 'D:/L1polymORF/Data/'
G1000SamplePath <- 'D:/L1polymORF/Data/1000GenomeSampleInfo.txt'
GapPath         <- 'D:/L1polymORF/Data/Gap_hg19.txt'
L1GRPath        <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
L1RefRangePath  <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
ChrLPath        <- 'D:/L1polymORF/Data/ChromLengthsHg19.Rdata'
DNAseDataPath   <- 'D:/L1polymORF/Data/DNAseInfo.RData'
InputPath       <- 'D:/L1polymORF/Data/SingletonAnalysis_unphased.RData'
L1RefPath       <- 'D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv'
CpGPath         <- 'D:/L1polymORF/Data/CpG_hg19.txt'
RecombDataPath  <- "D:/L1polymORF/Data/hg19RecombRate.txt"
ChromHHPath     <- "D:/L1polymORF/Data/EncodeBroadHMM/ChromHMMcombined.txt"
SummaryGRPath   <- "D:/L1polymORF/Data/SummaryGR"
OutputPathGR      <- 'D:/L1polymORF/Data/GRangesSummaries_1Mb.RData'

# Number of info columns in vcf file
NrInfoCols   <- 9

# False discovery rate for selected L1
FDR <- 0.1

# Specify range width for DNAse analysis
RangeWidth <- 10^6

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("Loading and processing data ...")

# Load previously generated objects
load(L1GRPath)
load(ChrLPath)
load(L1RefRangePath)
load(InputPath)

# Make genomic ranges for L1SingletonCoeffs
L1SingletonCoeffs$chromosome <- paste("chr", L1SingletonCoeffs$Chrom, sep = "")
L1SingletonCoeffs_GR <- makeGRangesFromDataFrame(L1SingletonCoeffs, 
                                                 seqnames.field = "chromosome",
                                                 start.field = "Pos",
                                                 end.field = "Pos")

# Read information about 1000 genome samples
SampleInfo  <- read.table(G1000SamplePath, header = T)
SampleMatch <- match(SampleColumns, SampleInfo$sample)
Pops        <- SampleInfo$super_pop[SampleMatch]
NrS         <- length(SampleColumns)

# Table for each L1 how often it occurs in each population
UniquePops <- unique(SampleInfo$super_pop)
PopFreq <- sapply(UniquePops, function(x){
  blnPop <- Pops == x
  rowSums(L1_1000G[,SampleColumns[blnPop]])
})
colnames(PopFreq) <- UniquePops

# Match coefficients to 1000 genome data
ChromPosCoeff     <- paste(L1SingletonCoeffs$chromosome, L1SingletonCoeffs$Pos)
ChromPos1000G     <- paste(L1_1000G_reduced$chromosome, L1_1000G_reduced$POS)
MatchCoeff1000G   <- match(ChromPosCoeff, ChromPos1000G)
L1SingletonCoeffs <- cbind(L1SingletonCoeffs, PopFreq[MatchCoeff1000G,])

# Define more genomic ranges
GeneGR <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
ExonGR <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
PromGR <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream = 10000)
CDSGR  <- cds(TxDb.Hsapiens.UCSC.hg19.knownGene)
IntronGRList   <- intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                      use.names = T)
FiveUTRGRList  <- fiveUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                       use.names = T)
ThreeUTRGRList <- threeUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                        use.names = T)
sum(width(GeneGR)/10^6) / sum(ChromLengthsHg19/10^6)


# Read DNAse hypersensitivity data (generated in script CombineEncodeInfo_DNAse)
load(DNAseDataPath)

# Read in table with regulatory elements
RegTable      <- read.table(ChromHHPath, header = T)
blnAllCellTypes <- RegTable$CellType == 
  "Gm12878,H1hesc,Hepg2,Hmec,Hsmm,Huvec,K562,Nhek,Nhlf"

# Get indices of different regulatory elements
idxEnhancer <- grep("Enhancer", RegTable$name)
idxTxn      <- grep("Txn", RegTable$name)
idxProm     <- grep("Promoter", RegTable$name)
idxHetero   <- grep("Heterochrom", RegTable$name)
idxRepr     <- union(grep("Insulator", RegTable$name),
                     grep("Repressed", RegTable$name))

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
EnhancerGR  <- makeGRangesFromDataFrame(RegTable[idxEnhancer,])
TxnGR       <- makeGRangesFromDataFrame(RegTable[idxTxn,])
PromGR      <- makeGRangesFromDataFrame(RegTable[idxProm, ])
ReprGR      <- makeGRangesFromDataFrame(RegTable[idxRepr, ])
rm(list = "RegTable")
gc()

# Read in table with L1HS from the referemce genome
L1RefTab <- read.csv(L1RefPath)
L1RefGR <- makeGRangesFromDataFrame(L1RefTab, seqnames.field = "genoName",
                                    start.field = "genoStart",
                                    end.field = "genoEnd")

# Read in CpG data and turn it into GRanges
CpGtable <- read.delim(CpGPath, header = T)
CpGGR    <- makeGRangesFromDataFrame(CpGtable, start.field = "chromStart",
                                     end.field = "chromEnd")

# Read in file and create GRanges
RecData <- read.delim(RecombDataPath)
Rec_GR  <- makeGRangesFromDataFrame(RecData)

# Read data on conserved 
PhastCons <- read.delim("D:/L1polymORF/Data/phastConsElements100way.txt", 
                        header = F)

# Add columns to 1000 genome data
L1_1000G$blnOLGene  <- overlapsAny(L1_1000G_GR_hg19, GeneGR)
L1_1000G$blnOLProm  <- overlapsAny(L1_1000G_GR_hg19, PromGR)
L1_1000G$blnOLExon  <- overlapsAny(L1_1000G_GR_hg19, ExonGR)
L1_1000G$blnOLCpG   <- overlapsAny(L1_1000G_GR_hg19, CpGGR)
L1_1000G$Dist2CpG   <- Dist2Closest(L1_1000G_GR_hg19, CpGGR)
L1_1000G$L1StartNum <- as.numeric(as.character(L1_1000G$L1Start))
L1_1000G$L1EndNum   <- as.numeric(as.character(L1_1000G$L1End))
L1_1000G$blnFull    <- L1_1000G$L1StartNum <= 1 & L1_1000G$L1EndNum >= 6000

cat("done!\n")


###################################################
#                                                 #
#     Summarize data per genomic range            #
#                                                 #
###################################################

cat("Summarizing data per genomic range ...")

# Summarize DNAse hypesensitivity data per RangeWidth
DNAse2Summarize  <- cbind(DNAseData[,c("chrom", "start", "end")], 
                          ScoreByStem)
DNAseSummaryList <- AggregateDataPerGRanges(
  Data2Summarize = DNAse2Summarize,
  Cols2Summarize = colnames(ScoreByStem),
  SeqNameCol = "chrom",
  StartCol = "start",
  EndCol = "end",
  RangeWidth = RangeWidth, 
  ChromLengths = ChromLengthsHg19,
  blnAddGRInfo = T,
  blnReturnDFOnly = F,
  NoValue = 0
)

# Create summary genomic ranges and data per summary range
SummaryGR        <- DNAseSummaryList$SummaryGR
DataPerSummaryGR <- DNAseSummaryList$SummarizedData
rm(list = "DNAseSummaryList")
gc()

###################
# Add more variables to summary
###################

# Add overlap counts
DataPerSummaryGR$L1Count <- countOverlaps(SummaryGR, L1_1000G_GR_hg19)
DataPerSummaryGR$L1Count_Full <- countOverlaps(SummaryGR, 
                                  L1_1000G_GR_hg19[which(L1_1000G$blnFull)])

# Add other L1 info to summary
L1Data2Summarize <- L1_1000G[ ,c("CHROM", "POS", "L1StartNum", "InsLength",
                                 "blnFull")]
L1Data2Summarize$CHROM <- paste("chr", L1Data2Summarize$CHROM, sep = "")
L1SummaryDf <- AggregateDataPerGRanges(
  Data2Summarize = L1Data2Summarize,
  SummaryGR = SummaryGR,
  Cols2Summarize = c("L1StartNum", "InsLength", "blnFull"),
  SeqNameCol = "CHROM",
  StartCol = "POS",
  EndCol = "POS",
  RangeWidth = RangeWidth, 
  ChromLengths = ChromLengthsHg19,
  blnAddGRInfo = T)

# Merge existing data with new ones
DataPerSummaryGR <- merge(DataPerSummaryGR, L1SummaryDf, sort = F)
rm(list = "L1SummaryDf")
gc()


# Add GC to summary
cat("calculating GC content ...")
SummaryViews <- BSgenomeViews(BSgenome.Hsapiens.UCSC.hg19, SummaryGR)
GCFreq       <- letterFrequency(SummaryViews, letters = "GC")
DataPerSummaryGR$GC <- GCFreq
cat("done\n")

# Add target site frequency to summary
cat("calculating target site and telomere hexamere frequency ...")
SeqWindows <- getSeq(BSgenome.Hsapiens.UCSC.hg19, SummaryGR)
TTTTAAFreq   <- vcountPattern("TTTTAA", SeqWindows)
TTAAAAFreq   <- vcountPattern("TTAAAA", SeqWindows)
DataPerSummaryGR$TargetFreq <- TTTTAAFreq + TTAAAAFreq

# Add telomere-containing hexamere frequency to summary
TTAGGGFreq    <- vcountPattern("TTAGGG", SeqWindows)
TTAGGGFreq_RC <- vcountPattern("CCCTAA", SeqWindows)
DataPerSummaryGR$TeloFreq <- TTAGGGFreq + TTAGGGFreq_RC
reverseComplement(DNAString('TTAGGG'))

cat("done\n")

cat("adding more columns (CpG, recombination frequency, phastCons) ...")
# Add CpG data
DataPerSummaryGR$CpGCount <- countOverlaps(SummaryGR, CpGGR)
CpGSummaryDF <- AggregateDataPerGRanges(
  Data2Summarize = CpGtable,
  Cols2Summarize = c("length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp"),
  SeqNameCol     = "chrom",
  StartCol       = "chromStart",
  EndCol         = "chromEnd",
  RangeWidth     = RangeWidth, 
  ChromLengths = ChromLengthsHg19,
  blnAddGRInfo = T)
CpGSummaryDF[is.na(CpGSummaryDF$length),c("length", "cpgNum", "gcNum")] <- 0

# Merge existing data with new ones
DataPerSummaryGR <- merge(DataPerSummaryGR, CpGSummaryDF, sort = F)
rm(list = c("CpGSummaryDF", "CpGGR", "CpGtable"))
gc()

# Add recombination data
RecSummaryDF <- AggregateDataPerGRanges(
  Data2Summarize = RecData,
  Cols2Summarize = c("decodeAvg",  "decodeFemale", "decodeMale", 
                     "marshfieldAvg", "marshfieldMale", "genethonAvg",
                     "genethonFemale", "genethonMale"),
  SeqNameCol     = "chrom",
  StartCol       = "chromStart",
  EndCol         = "chromEnd",
  RangeWidth     = RangeWidth, 
  ChromLengths = ChromLengthsHg19,
  blnAddGRInfo = T)

# Merge existing data with new ones
DataPerSummaryGR <- merge(DataPerSummaryGR, RecSummaryDF, sort = F)
rm(list = c("RecSummaryDF", "RecData"))
gc()

# Add conservation data data
PhastConsSummaryDF <- AggregateDataPerGRanges(
  Data2Summarize = PhastCons,
  Cols2Summarize = "V6",
  SeqNameCol     = "V2",
  StartCol       = "V3",
  EndCol         = "V4",
  RangeWidth     = RangeWidth, 
  ChromLengths = ChromLengthsHg19,
  blnAddGRInfo = T)

# Merge existing data with new ones
DataPerSummaryGR <- merge(DataPerSummaryGR, PhastConsSummaryDF, sort = F)
rm(list = c("PhastConsSummaryDF", "PhastCons"))
gc()
colnames(DataPerSummaryGR)[colnames(DataPerSummaryGR) == "V6"] <- "phastCons"

# Summarize singleton coefficient per 1 Mb window
CoeffSummaryDF <- AggregateDataPerGRanges(
  Data2Summarize = L1SingletonCoeffs,
  Cols2Summarize = c("coef",  "se.coef."),
  SeqNameCol     = "chromosome",
  StartCol       = "Pos",
  EndCol         = "Pos",
  RangeWidth     = RangeWidth, 
  ChromLengths = ChromLengthsHg19,
  blnAddGRInfo = T)

# Merge existing data with new ones
DataPerSummaryGR <- merge(DataPerSummaryGR, CoeffSummaryDF, sort = F)
rm(list = c("CoeffSummaryDF"))
gc()
colnames(DataPerSummaryGR)[colnames(DataPerSummaryGR) == "coef"] <- "SinglCoef"

# Add various counts (better: calculate sum of width)
CoverGR <- GeneGR
CalcGRCoverage <- function(SummaryGR, CoverGR){
  WidthDF <- data.frame(Chrom = as.vector(seqnames(CoverGR)),
    Start = start(CoverGR), End = end(CoverGR), Width = width(CoverGR))
  AggWidth <- AggregateDataPerGRanges(
    SummaryGR = SummaryGR,
    Data2Summarize = WidthDF,
    Cols2Summarize = c("Width"),
    SeqNameCol     = "Chrom",
    StartCol       = "Start",
    EndCol         = "End",
    ChromLengths = ChromLengthsHg19,
    SummaryFun = sum,
    blnAddGRInfo = T)
  AggWidth$Width
}
DataPerSummaryGR$EnhancerCount <- CalcGRCoverage(SummaryGR, EnhancerGR)
DataPerSummaryGR$TxnCount      <- CalcGRCoverage(SummaryGR, TxnGR)
DataPerSummaryGR$PromCount     <- CalcGRCoverage(SummaryGR, PromGR)
DataPerSummaryGR$ReprCount     <- CalcGRCoverage(SummaryGR, ReprGR)
DataPerSummaryGR$GeneCount     <- CalcGRCoverage(SummaryGR, GeneGR)

# Get a histogram of L1 frequency
hist(L1_1000G$Frequency, breaks = seq(0, 1, 10^(-3)), xlim = c(0, 0.1))

# Group L1 ranges in three frequency groups
L1_1000G_GR_hg19_low  <- L1_1000G_GR_hg19[L1_1000G$Frequency <= 10^(-3)]
L1_1000G_GR_hg19_med  <- L1_1000G_GR_hg19[L1_1000G$Frequency > 10^(-3) &
                                            L1_1000G$Frequency <= 0.02]
L1_1000G_GR_hg19_high <- L1_1000G_GR_hg19[L1_1000G$Frequency  > 0.02]

# Add overlap counts for different L1 frequencies
DataPerSummaryGR$L1Count_low  <- countOverlaps(SummaryGR, L1_1000G_GR_hg19_low)
DataPerSummaryGR$L1Count_med  <- countOverlaps(SummaryGR, L1_1000G_GR_hg19_med)
DataPerSummaryGR$L1Count_high <- countOverlaps(SummaryGR, L1_1000G_GR_hg19_high)

# Add overlap counts for references L1
DataPerSummaryGR$L1Count_ref  <- countOverlaps(SummaryGR, L1RefGR)

cat("done!\n")

###################################################
#                                                 #
#     Save data                                   #
#                                                 #
###################################################

cat("Saving data ...")
save(file = OutputPathGR, list = c("DataPerSummaryGR", "RangeWidth", "SummaryGR"))

# Export tables used for 1000 genome summaries
for (i in 1:4){
  End <- min(i*1000, length(SummaryGR))
  ExportGR <- SummaryGR[(1000*(i - 1) + 1):End]
  exportPath <- paste(SummaryGRPath, i, sep = "")
  ExportTab <- cbind(as.vector(seqnames(ExportGR)),
                     start(ExportGR), end(ExportGR))
  write.table(ExportTab, exportPath, quote = F, row.names = F,
              col.names = F)
  #export.bed(ExportGR, con = exportPath)
}

cat("done!\n")



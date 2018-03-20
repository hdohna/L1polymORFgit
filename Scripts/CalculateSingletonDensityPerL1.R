# The script below reads calculates the singleton density per L1
# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Specify file paths
DataPath <- 'D:/L1polymORF/Data/'
G1000SamplePath <- 'D:/L1polymORF/Data/1000GenomeSampleInfo.txt'
L1GRPath        <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
L1RefRangePath  <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'

# Specify parameters
NrInfoCols <- 9

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

# Load previously generated objects
load(L1GRPath)

##########
# 
##########

# Specify chromosome
Chr <- 9
cat("********   Analyzing chromosome,", Chr, "    **********\n")

# Read singleton file
cat("Reading singleton file ...")
SingletonPath <- paste(DataPath, "Singleton_SNP_chr", Chr, sep = "")
Singletons    <- read.table(SingletonPath)
SCols         <- GetSingletonColumns(Singletons)
SingletonGR   <- makeGRangesFromDataFrame(Singletons, seqnames.field = "V1",
                                          start.field = "V2",
                                          end.field = "V2")
cat("Done!\n")

# Read LINE-1 vcf file
cat("Reading LINE-1 vcf file ...")
Line1VcfPath <- paste(DataPath, "LINE1chr", Chr, ".vcf", sep = "")
Line1Vcf     <- read.table(Line1VcfPath)
cat("Done!\n")

# Subset LINE-1 vcf file
blnSampleCols <- colnames(L1_1000G) %in% SampleColumns
NrWith <- apply(Line1Vcf[,blnSampleCols], 1, function(x) length(grep("1", x)))
hist(NrWith, breaks = seq(0, 3000, 5), xlim = c(0, 300))
idx5   <- which(NrWith >= 5)
idxColHomoWith    <- which(L1_1000G[idxChr[i], SampleColumns] == 2)
idxColHHomoWithout <- which(L1_1000G[idxChr[i], SampleColumns] == 0)

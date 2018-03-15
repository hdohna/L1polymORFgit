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

# Load singleton file
SingletonPath <- paste("D:/L1polymORF/Data/Singleton_SNP_chr", Chr, sep = "")
Singletons    <- read.table(SingletonPath)
Singletons$Col <- ceiling((Singletons$V10 %% 10016)/4)

# Subset L1_1000G
idxChr <- which(L1_1000G$CHROM == Chr)
NrHomoZygWith   <- rowSums(L1_1000G[idxChr, SampleColumns] == 2)
idx5 <- which(NrHomoZygWith >= 5)
i <- idx5[1]
idxColHWith    <- which(L1_1000G[idxChr[i], SampleColumns] == 2)
idxColHWithout <- which(L1_1000G[idxChr[i], SampleColumns] == 0)

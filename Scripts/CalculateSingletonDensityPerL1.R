# The script below reads calculates the singleton density per L1
##########################################
#                                        #
#     Load packages                      #
#                                        #
##########################################

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(survival)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

# Specify file paths
DataPath <- 'D:/L1polymORF/Data/'
G1000SamplePath <- 'D:/L1polymORF/Data/1000GenomeSampleInfo.txt'
L1GRPath        <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
L1RefRangePath  <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
ChrLPath        <- 'D:/L1polymORF/Data/ChromLengthsHg19.Rdata'

# Number of info columns in vcf file
NrInfoCols <- 9

# Minimum number of carriers for a LINE-1 to be analyzed
MinNrCarrier <- 5

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

# Load previously generated objects
load(L1GRPath)
load(ChrLPath)

##########
# 
##########

# Specify chromosome
Chr <- 9
cat("********   Analyzing chromosome,", Chr, "    **********\n")

# Get chromosome length
ChrL <- ChromLengthsHg19[Chr]

# Read singleton file
cat("Reading singleton file ...")
SingletonPath  <- paste(DataPath, "Singleton_SNP_chr", Chr, sep = "")
Singletons     <- read.table(SingletonPath)
SCols          <- GetSingletonColumns(Singletons)
SCols_rev      <- SCols[nrow(SCols):1,]
Singletons_rev <- Singletons[nrow(Singletons):1, ]
cat("Done!\n")

# Read LINE-1 vcf file
cat("Reading LINE-1 vcf file ...")
Line1VcfPath <- paste(DataPath, "LINE1chr", Chr, ".vcf", sep = "")
Line1Vcf     <- read.table(Line1VcfPath, as.is = T)
cat("Done!\n")

# Subset LINE-1 vcf file
blnSampleCols <- colnames(L1_1000G) %in% SampleColumns
idxSampleCols <- which(blnSampleCols)

# Get for each L1 the number of carriers
NrCarriers  <- apply(Line1Vcf[,blnSampleCols], 1, function(x) length(grep("1", x)))
idxEnough   <- which(NrCarriers >= MinNrCarrier)

# Initialize objects
blnSingl1 <- which(SCols$Allele == 1) 
blnSingl3 <- which(SCols$Allele == 3)
idxSCols  <- 1:length(SampleColumns)
L1Col     <- rep(0, length(SampleColumns))

# Loop over L1 with enough carriers and estimate effect of L1 on singleton
# densities
L1Coeffs <- sapply (idxEnough, function(i) {
  cat("Processing L1", which(idxEnough == i), "out of", length(idxEnough), "\n")
  blnBeforeL1 <- Singletons_rev$V2 < Line1Vcf$V2[i]
  blnAfterL1  <- Singletons$V2 > Line1Vcf$V2[i]
  idxCarrier1 <- grep("1\\|", Line1Vcf[i,blnSampleCols])
  idxCarrier3 <- grep("\\|1", Line1Vcf[i,blnSampleCols])
  CarrierMatch1After  <- match(idxSCols, SCols$Col[blnSingl1 & blnAfterL1])
  CarrierMatch1Before <- match(idxSCols, SCols_rev$Col[blnSingl1 & blnBeforeL1])
  CarrierMatch3After  <- match(idxSCols, SCols$Col[blnSingl3 & blnAfterL1])
  CarrierMatch3Before <- match(idxSCols, SCols_rev$Col[blnSingl3 & blnBeforeL1])
  DistBefore1 <- Line1Vcf$V2[i] - Singletons_rev$V2[CarrierMatch1Before]
  DistAfter1  <- Singletons_rev$V2[CarrierMatch1Before] - Line1Vcf$V2[i]
  DistBefore3 <- Line1Vcf$V2[i] - Singletons$V2[CarrierMatch3Before]
  DistAfter3  <- Singletons$V2[CarrierMatch3Before] - Line1Vcf$V2[i]

  L1_1 <- rep(0, length(SampleColumns))
  L1_1[idxCarrier1] <- 1
  L1_3 <- rep(0, length(SampleColumns))
  L1_3[idxCarrier3] <- 1
  CensorBefore <- rep(1, length(SampleColumns))
  CensorAfter  <- rep(1, length(SampleColumns))
  
  DistBefore  <- c(DistBefore1, DistBefore3)
  blnNABefore <- is.na(DistBefore)
  CensorBefore[blnNABefore] <- 0
  DistBefore[blnNABefore]   <- Line1Vcf$V2[i]
  
  DistAfter  <- c(DistAfter1, DistAfter3)
  blnNAAfter <- is.na(DistAfter)
  CensorAfter[blnNAAfter] <- 0
  DistAfter[blnNAAfter]   <- ChrL - Line1Vcf$V2[i]
  
  SurvObj <- Surv(time =  c(DistBefore, DistAfter),
                  event = c(CensorBefore, CensorAfter))
  L1         <- c(L1_1, L1_3, L1_1, L1_3)
  Direction  <- c(rep("Before", length(DistBefore)), 
                  rep("After",  length(DistAfter)))
  SampleIDs  <- rep(idxSCols, 4)

  CPH <- coxph(SurvObj ~ L1 + strata(Direction) + cluster(SampleIDs))
  SU <- summary(CPH)
  SU$coefficients
}) 

hist(L1Coeffs[1,])
mean(L1Coeffs[1,])
sqrt(var(L1Coeffs[1,]) / ncol(L1Coeffs))

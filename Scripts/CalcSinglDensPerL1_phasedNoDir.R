# The script below calculates the singleton density per L1 using phasing 
# information but not distinguishing between directions along the genome

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
DataPath        <- 'D:/L1polymORF/Data/'
G1000SamplePath <- 'D:/L1polymORF/Data/1000GenomeSampleInfo.txt'
GapPath         <- 'D:/L1polymORF/Data/Gap_hg19.txt'
L1GRPath        <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
L1RefRangePath  <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
ChrLPath        <- 'D:/L1polymORF/Data/ChromLengthsHg19.Rdata'
HiCFolderPath   <- 'D:/L1polymORF/Data/HiCData/'
OutputPath      <- 'D:/L1polymORF/Data/SingletonAnalysis_phasedNoDir.RData'

# Number of info columns in vcf file
NrInfoCols   <- 9

# Minimum number of carriers for a LINE-1 to be analyzed
MinNrCarrier <- 3

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

# Load previously generated objects
load(L1GRPath)
load(ChrLPath)

# Read information about 1000 genome samples
SampleInfo  <- read.table(G1000SamplePath, header = T)
SampleMatch <- match(SampleColumns, SampleInfo$sample)
Pops        <- SampleInfo$super_pop[SampleMatch]
NrS         <- length(SampleColumns)

##########
# Loop over chromosomes and estimate coefficients
##########

# Initialize data.frame for coefficients associated with L1
L1SingletonCoeffs <- data.frame()
NrS               <- length(SampleColumns)
idxSCols          <- 1:NrS 

# Specify chromosome
for (Chr in c(1:22, "X")) {

  cat("********   Analyzing chromosome,", Chr, "    **********\n")
  # Get chromosome length
  ChrL <- ChromLengthsHg19[Chr]
  
  # Read singleton file
  cat("Reading singleton file ...")
  SingletonPath  <- paste(DataPath, "Singleton_SNP_chr", Chr, sep = "")
  Singletons     <- read.table(SingletonPath)
  SCols          <- GetSingletonColumns(Singletons)
  SCols$Allele[SCols$Allele == 0] <- 1
  SCols$Allele[SCols$Allele == 2] <- 3
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
  blnSingl1 <- SCols$Allele == 1
  blnSingl3 <- SCols$Allele == 3
  
  # Loop over L1 with enough carriers and estimate effect of L1 on singleton
  # densities
  i <- idxEnough[1]
  L1Coeffs <- sapply (idxEnough, function(i) {
    cat("Processing L1", which(idxEnough == i), "out of", length(idxEnough), "\n")
    
    # Indicators for singletons that are before and after current L1
    blnBeforeL1 <- Singletons_rev$V2 < Line1Vcf$V2[i]
    blnAfterL1  <- Singletons$V2     > Line1Vcf$V2[i]
    
    # Index for L1 carriers on first abd second allele
    idxCarrier1 <- grep("1\\|", Line1Vcf[i,blnSampleCols])
    idxCarrier3 <- grep("\\|1", Line1Vcf[i,blnSampleCols])
    
    # Match singleton entries for all samples
    idxA1 <- which(blnSingl1 & blnAfterL1)
    idxB1 <- which(blnSingl1 & blnBeforeL1)
    idxA3 <- which(blnSingl3 & blnAfterL1)
    idxB3 <- which(blnSingl3 & blnBeforeL1)
    SampleMatch1After  <- match(idxSCols, SCols$Col[idxA1])
    SampleMatch1Before <- match(idxSCols, SCols_rev$Col[idxB1])
    SampleMatch3After  <- match(idxSCols, SCols$Col[idxA3])
    SampleMatch3Before <- match(idxSCols, SCols_rev$Col[idxB3])
    
    # Determine distance from L1 to singleton
    DistBefore1 <- Line1Vcf$V2[i] - Singletons_rev$V2[idxB1[SampleMatch1Before]]
    DistAfter1  <- Singletons$V2[idxA1[SampleMatch1After]] - Line1Vcf$V2[i]
    DistBefore3 <- Line1Vcf$V2[i] - Singletons_rev$V2[idxB3[SampleMatch3Before]]
    DistAfter3  <- Singletons$V2[idxA3[SampleMatch3After]] - Line1Vcf$V2[i]
    
    # Create indicator for L1 presence
    L1_1 <- rep(0, NrS)
    L1_1[idxCarrier1] <- 1
    L1_3 <- rep(0, NrS)
    L1_3[idxCarrier3] <- 1
    
    # Create censoring column
    Censor1 <- rep(1, length(SampleColumns))
    Censor3 <- rep(1, length(SampleColumns))
    
    # Replace NA distances by distance to either end of chromosome
    blnNABefore1 <- is.na(DistBefore1)
    blnNABefore3 <- is.na(DistBefore3)
    blnNAAfter1  <- is.na(DistAfter1)
    blnNAAfter3  <- is.na(DistAfter3)
    Censor1[blnNABefore1 | blnNAAfter1]    <- 0
    Censor3[blnNABefore3 | blnNAAfter3]    <- 0

    DistBefore1[blnNABefore1] <- Line1Vcf$V2[i]
    DistBefore3[blnNABefore3] <- Line1Vcf$V2[i]
    DistAfter1[blnNAAfter1]   <- ChrL - Line1Vcf$V2[i]
    DistAfter3[blnNAAfter3]   <- ChrL - Line1Vcf$V2[i]

    # Create two distances for the two alleles
    Dist1 <- DistBefore1 + DistAfter1
    Dist3 <- DistBefore3 + DistAfter3
    
    # Create survival object and perform Cox proportional hazards
    SurvObj <- Surv(time =  c(Dist1, Dist3), event = c(Censor1, Censor3))
    L1      <- c(L1_1, L1_3)
    CPH     <- coxph(SurvObj ~ L1 + rep(Pops, 2))
    SU      <- summary(CPH)
    Coeffs  <- SU$coefficients["L1",]
    names(Coeffs) <- colnames(SU$coefficients)
    Coeffs
  }) 
  
  # Create a data frame with L1 coefficients
  L1Coeffs <- data.frame(t(L1Coeffs))
  L1Coeffs$Chrom <- Chr
  L1Coeffs$Pos   <- Line1Vcf$V2[idxEnough]
  
  CHrPos_all   <- paste(L1_1000G$CHROM, L1_1000G$POS)
  CHrPos_match <- match(paste(Chr, Line1Vcf$V2[idxEnough]), CHrPos_all)
  L1Coeffs$InsLength <- L1_1000G_reduced$InsLength[CHrPos_match]
  L1Coeffs$Freq      <- L1_1000G_reduced$Frequency[CHrPos_match]
  L1Coeffs$L1Start   <- L1_1000G_reduced$L1Start[CHrPos_match]
  L1Coeffs$L1End     <- L1_1000G_reduced$L1End[CHrPos_match]
  L1Coeffs$L1Strand  <- L1_1000G_reduced$L1Strand[CHrPos_match]
  
  L1SingletonCoeffs <- rbind(L1SingletonCoeffs, L1Coeffs)

}

##########################
#                        #
#    Save image          #
#                        #
##########################

save.image(OutputPath)

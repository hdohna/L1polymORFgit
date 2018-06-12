# The script below reads calculates the singleton density per L1 by simulating
# low frequency L1

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

# Number of info columns in vcf file
NrInfoCols   <- 9

# Minimum number of carriers for a LINE-1 to be analyzed
MinNrCarrier  <- 10
NrCarrier2Sim <- 3

# Sample factor per chromosome
SampleFact <- 5

# Create an output path
OutputPath  <- paste('D:/L1polymORF/Data/SingletonAnalysis_SampledL1_Freq',
                     NrCarrier2Sim, '.RData', sep = "")


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
NrS <- length(SampleColumns)
idxSCols  <- 1:NrS 
Chr <- 1
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
  cat(length(idxEnough), "L1s occur in at least", MinNrCarrier, "carriers\n")
  idxSampledL1 <- sample(idxEnough, length(idxEnough) * SampleFact, replace = T)

  # Initialize objects
  blnSingl1 <- SCols$Allele == 1
  blnSingl3 <- SCols$Allele == 3
  
  # Loop over L1 with enough carriers and estimate effect of L1 on singleton
  # densities
  j <- 1
  idxSampledL1[1:10]
  L1Coeffs <- sapply (1:length(idxSampledL1), function(j) {
    cat("Processing L1", j, "out of", length(idxSampledL1), "\n")
    
    # Replace carriers by zero to generate low frequency L1
    i <- idxSampledL1[j]
    idxWith    <- grep("1",Line1Vcf[i, blnSampleCols])
    idxReplace <- sample(idxWith, length(idxWith) - NrCarrier2Sim)
    length(idxWith) - length(idxReplace)
    Line1Vcf[i, idxSampleCols[idxReplace]] <- '0|0'
    
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
    
    # Create an Allele column
    Alleles <- rep(c(rep(1, NrS ), rep(3, NrS )), 2)
    
    # Create censoring column
    CensorBefore <- rep(1, length(SampleColumns))
    CensorAfter  <- rep(1, length(SampleColumns))
    
    # Replace NA distances by distance to end
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
    SampleIDs  <- paste(rep(idxSCols, 4), Alleles, sep = "_")
    
    CPH <- coxph(SurvObj ~ L1 + strata(Direction) + cluster(SampleIDs) + rep(Pops, 4))
    SU <- summary(CPH)
    Coeffs <- SU$coefficients["L1",]
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
  L1Coeffs$Freq      <- NrCarrier2Sim / length(idxSampleCols)
  L1Coeffs$L1Start   <- L1_1000G_reduced$L1Start[CHrPos_match]
  L1Coeffs$L1End     <- L1_1000G_reduced$L1End[CHrPos_match]
  L1Coeffs$L1Strand  <- L1_1000G_reduced$L1Strand[CHrPos_match]
  
  L1SingletonCoeffs <- rbind(L1SingletonCoeffs, L1Coeffs)
#  plot(L1Coeffs$InsLength, L1Coeffs$coef, main = paste("chr", Chr))
  
}

# Assign table with L1 singleton coefficients to a new name and save
cat("Saving results into", OutputPath)
L1CoefName <- paste('L1SingletonCoeffs', NrCarrier2Sim, sep = "_")
assign(L1CoefName, L1SingletonCoeffs, envir = .GlobalEnv)
save(list = L1CoefName, file = OutputPath)

mean(L1SingletonCoeffs_3$Pr...z.., na.rm = T)

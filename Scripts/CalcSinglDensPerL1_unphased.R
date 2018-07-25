# The script below reads calculates the singleton density per L1 not using 
# phasing information

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
OutputPath      <- 'D:/L1polymORF/Data/SingletonAnalysis_unphased.RData'

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
Chr <- 1
# Specify chromosome
for (Chr in c(1:22, "X")) {
  cat("********   Analyzing chromosome", Chr, "    **********\n")
  
  # Get chromosome length and index of L1s for current chromosome
  ChrL   <- ChromLengthsHg19[Chr]
  idxChr <- which(L1_1000G$CHROM == Chr)
  
  # Read singleton file
  cat("Reading singleton file ...")
  # SingletonPath  <- paste(DataPath, "Singleton_SNP_chr", Chr, sep = "")
  # Singletons     <- read.table(SingletonPath)
  # SCols          <- GetSingletonColumns(Singletons)
  # SCols_rev      <- SCols[nrow(SCols):1,]
  # Singletons_rev <- Singletons[nrow(Singletons):1, ]

  SingletonPath  <- paste("D:/L1polymORF/Data/chr", Chr, ".singletons", sep = "")
  Singletons     <- read.table(SingletonPath, header = T)
  Singletons_rev <- Singletons[nrow(Singletons):1, ]
  
  cat("Done!\n")
  
  # Subset LINE-1 vcf file
  blnSampleCols <- colnames(L1_1000G) %in% SampleColumns
  idxSampleCols <- which(blnSampleCols)
  
  # Get for each L1 the number of carriers
  NrCarriers <- rowSums(L1_1000G[idxChr,idxSampleCols] > 0)
  idxEnough  <- idxChr[NrCarriers >= MinNrCarrier]
  
  # Loop over L1 with enough carriers and estimate effect of L1 on singleton
  # densities
  L1Coeffs <- sapply (idxEnough, function(i) {
    cat("Processing L1", which(idxEnough == i), "out of", length(idxEnough), "\n")
    
    # Indicators for singletons that are before and after current L1
    blnBeforeL1 <- Singletons_rev$POS < L1_1000G$POS[i]
    blnAfterL1  <- Singletons$POS     > L1_1000G$POS[i]
    
    # Index for L1 carriers on first and second allele
    idxCarrier  <- which(L1_1000G[i,blnSampleCols] > 0)

    # Match singleton entries for all samples
    idxA <- which(blnAfterL1)
    idxB <- which(blnBeforeL1)
    SampleMatchAfter  <- match(SampleColumns, Singletons$INDV[blnAfterL1])
    SampleMatchBefore <- match(SampleColumns, Singletons_rev$INDV[blnBeforeL1])

    # Determine distance from L1 to singleton
    DistBefore <- L1_1000G$POS[i] - Singletons_rev$POS[idxB[SampleMatchBefore]]
    DistAfter  <- Singletons$POS[idxA[SampleMatchAfter]] - L1_1000G$POS[i]

    # Create censoring column
    Censor <- rep(1, length(SampleColumns))

    # Replace NA distances by distance to end and set censor variable to zero
    blnNABefore <- is.na(DistBefore)
    blnNAAfter  <- is.na(DistAfter)
    DistBefore[blnNABefore] <- L1_1000G$POS[i]
    DistAfter[blnNAAfter]   <- ChrL - L1_1000G$POS[i]
    Censor[blnNAAfter | blnNABefore] <- 0
    
    # Add upstream and downstream distance to next singleton
    Dist <- DistBefore + DistAfter

    # Estimate cox ph model to distance to singletons
    SurvObj <- Surv(time = Dist,event = Censor)
    L1      <- t(L1_1000G[i,SampleColumns])
    CPH     <- coxph(SurvObj ~ L1 + Pops)
    SU      <- summary(CPH)
    Coeffs  <- SU$coefficients["L1",]
    names(Coeffs) <- colnames(SU$coefficients)
    Coeffs
  }) 
  
  # Create a data frame with L1 coefficients
  L1Coeffs           <- data.frame(t(L1Coeffs))
  L1Coeffs$Chrom     <- Chr
  L1Coeffs$Pos       <- L1_1000G$POS[idxEnough]
  L1Coeffs$InsLength <- L1_1000G_reduced$InsLength[idxEnough]
  L1Coeffs$Freq      <- L1_1000G_reduced$Frequency[idxEnough]
  L1Coeffs$L1Start   <- L1_1000G_reduced$L1Start[idxEnough]
  L1Coeffs$L1End     <- L1_1000G_reduced$L1End[idxEnough]
  L1Coeffs$L1Strand  <- L1_1000G_reduced$L1Strand[idxEnough]
  
  L1SingletonCoeffs <- rbind(L1SingletonCoeffs, L1Coeffs)
}

##########################
#                        #
#    Save image          #
#                        #
##########################

cat("*********  Saving results ... ")
save.image(OutputPath)
cat("done!  *********")

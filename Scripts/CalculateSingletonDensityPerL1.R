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
GapPath         <- 'D:/L1polymORF/Data/Gap_hg19.txt'
L1GRPath        <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
L1RefRangePath  <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
ChrLPath        <- 'D:/L1polymORF/Data/ChromLengthsHg19.Rdata'
HiCFolderPath   <- 'D:/L1polymORF/Data/HiCData/'
OutputPath     <- 'D:/L1polymORF/Data/SingletonAnalysis.RData'

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
NrS <- length(SampleColumns)
idxSCols  <- 1:NrS 

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
  L1Coeffs$Freq      <- L1_1000G_reduced$Frequency[CHrPos_match]
  
  L1SingletonCoeffs <- rbind(L1SingletonCoeffs, L1Coeffs)
#  plot(L1Coeffs$InsLength, L1Coeffs$coef, main = paste("chr", Chr))
  
}

##########################################
#                                        #
#        Add columns                     #
#                                        #
##########################################

# Indicator for full-length
L1SingletonCoeffs$blnFull <- L1SingletonCoeffs$InsLength >= 6000

# Indicator for significant effect
L1SingletonCoeffs$blnSig <- L1SingletonCoeffs$Pr...z.. < 0.05

# Indicator for  selection
L1SingletonCoeffs$blnSelect <- L1SingletonCoeffs$Pr...z.. < 0.05 &
  L1SingletonCoeffs$coef < 0

# Bin for insertion length
L1SingletonCoeffs$InsLBins <- cut(L1SingletonCoeffs$InsLength, 
                                  breaks = seq(0, 6500, 500))

##########################################
#                                        #
#        Plot coefficients               #
#                                        #
##########################################

# Plot coefficients vs insertion length
plot(L1SingletonCoeffs$InsLength, L1SingletonCoeffs$coef, 
     xlab = "Insertion length", ylab = "Singleton coefficient")
with(L1SingletonCoeffs, points(InsLength[blnSelect], coef[blnSelect], 
       col = "red"))
MeanCoeffPerL <- aggregate(L1SingletonCoeffs[,c("coef", "InsLength", "blnSelect")],
                           by = list(L1SingletonCoeffs$InsLBins), 
                           FUN = function(x) mean(x, na.rm = T))
plot(MeanCoeffPerL$InsLength, MeanCoeffPerL$coef)
plot(MeanCoeffPerL$InsLength, MeanCoeffPerL$blnSelect)

##########################################
#                                        #
#        Regress coefficients            #
#                                        #
##########################################

# Make genomic ranges for L1SingletonCoeffs
L1SingletonCoeffs$chromosome <- paste("chr", L1SingletonCoeffs$Chrom, sep = "")
L1SingletonCoeffs_GR <- makeGRangesFromDataFrame(L1SingletonCoeffs, 
                                                 seqnames.field = "chromosome",
                                                 start.field = "Pos",
                                                 end.field = "Pos")

# Caclulate distance to genes
L1SingletonCoeffs$Dist2Gene <- Dist2Closest(L1SingletonCoeffs_GR, 
                                            genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
#######
# Regress against length
#######

LM_All_Interact_binom <- glm(1*blnSelect ~ InsLength + blnFull + InsLength*blnFull, 
                             data = L1SingletonCoeffs, 
                             subset =robust.se > 0 & (!is.na(blnSelect)),
                             weights = 1/L1SingletonCoeffs$robust.se, family  = binomial)
summary(LM_All_Interact_binom)
table(1*L1SingletonCoeffs$blnSelect)

PropSmoothed_InsL <- supsmu(L1SingletonCoeffs$InsLength,
                            1*L1SingletonCoeffs$blnSelect,
                            wt = 1/L1SingletonCoeffs$robust.se)
plot(PropSmoothed_InsL$x, PropSmoothed_InsL$y, xlab = "Insertion length [bp]",
     ylab = "Proportion of L1 with selection signal", type = "l",
     ylim = c(0, 0.2), col = "red")
lines(PropSmoothed_InsL$x[PropSmoothed_InsL$x>=6000], 
      PropSmoothed_InsL$y[PropSmoothed_InsL$x>=6000])
legend("topright", legend = c("Full-length L1", "Fragment L1"), col = c("black", "red"), 
       lty = c(1, 1))
CreateDisplayPdf('D:/L1polymORF/Figures/PropSelectVsInsLength.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

#######
# Regress against distance to genes
#######

# Regress coefficient against distance to genes
LM_Full <- lm(coef ~ Dist2Gene, data = L1SingletonCoeffs, 
         subset = L1SingletonCoeffs$InsLength >= 6000 & L1SingletonCoeffs$robust.se > 0,
         weights = 1/L1SingletonCoeffs$robust.se)
summary(LM_Full)

LM_Fragm <- lm(coef ~ Dist2Gene, data = L1SingletonCoeffs, 
              subset = L1SingletonCoeffs$InsLength < 5900 & L1SingletonCoeffs$robust.se > 0,
              weights = 1/L1SingletonCoeffs$robust.se)
summary(LM_Fragm)

LM_All <- lm(coef ~ Dist2Gene, data = L1SingletonCoeffs, 
              subset = L1SingletonCoeffs$robust.se > 0,
              weights = 1/L1SingletonCoeffs$robust.se)
LM_All <- lm(coef ~ Dist2Gene, data = L1SingletonCoeffs)
summary(LM_All)
sum(L1SingletonCoeffs$robust.se == 0)

LM_All_Interact_binom <- glm(1*blnSelect ~ Dist2Gene + blnFull + blnFull*Dist2Gene, 
                             data = L1SingletonCoeffs, 
                             subset =robust.se > 0 & (!is.na(blnSelect)),
                             weights = 1/L1SingletonCoeffs$robust.se, family  = binomial)
summary(LM_All_Interact_binom)


##########
# Process loops
##########

# Auxiliary unction to create genomic ranges for both sides of a loop
getLoopGRs <- function(FileName){
  Loops <- read.delim(FileName)
  Loops$chr1 <- paste("chr", Loops$chr1, sep = "")
  Loops$chr2 <- paste("chr", Loops$chr2, sep = "")
  LoopsGR1 <- makeGRangesFromDataFrame(Loops, seqnames.field = "chr1", 
                                       start.field="x1", end.field = "x2")
  LoopsGR2 <- makeGRangesFromDataFrame(Loops, seqnames.field = "chr2", 
                                       start.field="y1", end.field = "y2")
  blnOverlapLoop <- overlapsAny(LoopsGR1, LoopsGR2)
  AllLoops <- c(LoopsGR1, LoopsGR2[!blnOverlapLoop])
  GRangesList(LoopsGR1 = LoopsGR1, LoopsGR2 = LoopsGR2, AllLoops = AllLoops)
}

# List of files with domains and loops
DomainFiles <- list.files(HiCFolderPath, pattern = "domainlist.txt",
                          full.names = T)
DomainFiles <- DomainFiles[! DomainFiles %in% grep(".gz", DomainFiles, value = T)]
LoopFiles <- list.files(HiCFolderPath, pattern = "looplist.txt",
                        full.names = T)
LoopFiles <- LoopFiles[! LoopFiles %in% grep(".gz", LoopFiles, value = T)]

# Get list of domain and loop ranges
DomainGRList   <- lapply(DomainFiles, getLoopGRs)
LoopGRList     <- lapply(LoopFiles, getLoopGRs)
names(DomainGRList) <- sapply(DomainFiles, 
                              function(x) strsplit(x, "_")[[1]][2])
names(LoopGRList) <- sapply(LoopFiles, 
                            function(x) strsplit(x, "_")[[1]][2])

# Summarize domain ranges (union or intersect)
DomainGR_Intersect <- DomainGRList[[1]]$AllLoops
LoopGR_Intersect   <- LoopGRList[[1]]$AllLoops
for (i in 2:length(LoopGRList)){
  DomainGR_Intersect <- intersect(DomainGR_Intersect, DomainGRList[[i]]$AllLoops)
  LoopGR_Intersect   <- intersect(LoopGR_Intersect, LoopGRList[[i]]$AllLoops)
}

# Determine distance to loops and domains
L1SingletonCoeffs$Dist2Loop <- Dist2Closest(L1SingletonCoeffs_GR, 
                                            LoopGR_Intersect)
L1SingletonCoeffs$Dist2Domain <- Dist2Closest(L1SingletonCoeffs_GR, 
                                              DomainGR_Intersect)
L1SingletonCoeffs$Dist2LoopBins <- cut(L1SingletonCoeffs$Dist2Loop,
                                       breaks = seq(0, 24*10^6, 10^5))
max(L1SingletonCoeffs$Dist2Loop)

#######
# Regress against distance to loops
#######

# Regress coefficient against distance to loops 
LM_Full <- lm(coef ~ Dist2Loop, data = L1SingletonCoeffs, 
              subset = L1SingletonCoeffs$InsLength >= 6000 & L1SingletonCoeffs$robust.se > 0,
              weights = 1/L1SingletonCoeffs$robust.se)
summary(LM_Full)

LM_Fragm <- lm(coef ~ Dist2Loop, data = L1SingletonCoeffs, 
               subset = L1SingletonCoeffs$InsLength < 5900 & L1SingletonCoeffs$robust.se > 0,
               weights = 1/L1SingletonCoeffs$robust.se)
summary(LM_Fragm)

LM_All <- lm(coef ~ Dist2Loop, data = L1SingletonCoeffs, 
               subset = L1SingletonCoeffs$robust.se > 0,
               weights = 1/L1SingletonCoeffs$robust.se)
summary(LM_All)

LM_All_Interact <- lm(coef ~ Dist2Loop + blnFull + blnFull*Dist2Loop, data = L1SingletonCoeffs, 
             subset = L1SingletonCoeffs$robust.se > 0,
             weights = 1/L1SingletonCoeffs$robust.se)
summary(LM_All_Interact)

LM_All_Interact_binom <- glm(1*blnSelect ~ Dist2Loop + blnFull + blnFull*Dist2Loop, 
                             data = L1SingletonCoeffs, 
                      subset =robust.se > 0 & (!is.na(blnSelect)),
                      weights = 1/L1SingletonCoeffs$robust.se, family  = binomial)
summary(LM_All_Interact_binom)
table(1*L1SingletonCoeffs$blnSelect)

# Proportion of selected L1 per distance class
MeanSelectPerDist <- aggregate(L1SingletonCoeffs[,c("Dist2Loop", "blnSelect")],
                           by = list(L1SingletonCoeffs$Dist2LoopBins,
                                     L1SingletonCoeffs$blnFull), 
                           FUN = function(x) mean(x, na.rm = T))
with(MeanSelectPerDist, {
  plot(Dist2Loop[Group.2], blnSelect[Group.2]);
  points(Dist2Loop[!Group.2], blnSelect[!Group.2], col = "red");
  })
PropSmoothed_Full <- supsmu(L1SingletonCoeffs$Dist2Loop[L1SingletonCoeffs$blnFull],
                            1*L1SingletonCoeffs$blnSelect[L1SingletonCoeffs$blnFull],
                            wt = 1/L1SingletonCoeffs$robust.se[L1SingletonCoeffs$blnFull])
PropSmoothed_Fragm <- supsmu(L1SingletonCoeffs$Dist2Loop[!L1SingletonCoeffs$blnFull],
                            1*L1SingletonCoeffs$blnSelect[!L1SingletonCoeffs$blnFull],
                            wt = 1/L1SingletonCoeffs$robust.se[!L1SingletonCoeffs$blnFull])
plot(PropSmoothed_Full$x, PropSmoothed_Full$y, xlab = "Distance to closest loop [bp]",
     ylab = "Proportion of L1 with selection signal", type = "l",
     ylim = c(0, 0.2))
lines(PropSmoothed_Fragm$x, PropSmoothed_Fragm$y, col = "red")
# points(L1SingletonCoeffs$Dist2Loop[L1SingletonCoeffs$blnFull], 
#        L1SingletonCoeffs$blnSelect[L1SingletonCoeffs$blnFull])
# points(L1SingletonCoeffs$Dist2Loop[!L1SingletonCoeffs$blnFull], 
#        L1SingletonCoeffs$blnSelect[!L1SingletonCoeffs$blnFull], col = "red")
legend("topright", legend = c("Full-length L1", "Fragment L1"), col = c("black", "red"), 
       lty = c(1, 1))
CreateDisplayPdf('D:/L1polymORF/Figures/PropSelectVsDist2Loop.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

#######
# Regress against distance to loops and genes
#######

# Regress coefficient against distance to loops 
LM_Full <- lm(coef ~ Dist2Loop + Dist2Gene, data = L1SingletonCoeffs, 
              subset = L1SingletonCoeffs$InsLength >= 6000 & L1SingletonCoeffs$robust.se > 0,
              weights = 1/L1SingletonCoeffs$robust.se)
summary(LM_Full)

LM_Fragm <- lm(coef ~ Dist2Loop + Dist2Gene, data = L1SingletonCoeffs, 
               subset = L1SingletonCoeffs$InsLength < 5900 & L1SingletonCoeffs$robust.se > 0,
               weights = 1/L1SingletonCoeffs$robust.se)
summary(LM_Fragm)

##########################
#                        #
#    Save image          #
#                        #
##########################


save.image(OutputPath)

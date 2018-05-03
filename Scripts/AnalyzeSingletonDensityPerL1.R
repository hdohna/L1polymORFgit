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
DataPath        <- 'D:/L1polymORF/Data/'
G1000SamplePath <- 'D:/L1polymORF/Data/1000GenomeSampleInfo.txt'
GapPath         <- 'D:/L1polymORF/Data/Gap_hg19.txt'
L1GRPath        <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
L1RefRangePath  <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
ChrLPath        <- 'D:/L1polymORF/Data/ChromLengthsHg19.Rdata'
HiCFolderPath   <- 'D:/L1polymORF/Data/HiCData/'
InputPath       <- 'D:/L1polymORF/Data/SingletonAnalysis.RData'

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
load(InputPath)

# Read information about 1000 genome samples
SampleInfo  <- read.table(G1000SamplePath, header = T)
SampleMatch <- match(SampleColumns, SampleInfo$sample)
Pops        <- SampleInfo$super_pop[SampleMatch]
NrS         <- length(SampleColumns)

##########################################
#                                        #
#        Add columns                     #
#                                        #
##########################################

# Indicator for full-length
L1SingletonCoeffs$blnFull <- L1SingletonCoeffs$InsLength >= 6000
sum(L1SingletonCoeffs$InsLength <= 100)
# Indicator for significant effect
L1SingletonCoeffs$blnSig <- p.adjust(L1SingletonCoeffs$Pr...z..) < 0.1

hist(L1SingletonCoeffs$Pr...z..[L1SingletonCoeffs$coef < 0], 
     breaks = seq(0, 1, 0.01))

# Indicator for  selection
L1SingletonCoeffs$blnSelect <- L1SingletonCoeffs$blnSig &
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

# Subset L1 coefficients
sum(is.na(L1SingletonCoeffs$InsLength))
L1SingletonCoeffs_subset <- subset(L1SingletonCoeffs, 
                                   subset = robust.se > 0 & (!is.na(blnSelect)))
LM_All_Interact_binom <- glm(1L*blnSelect ~ InsLength + blnFull, 
                             data = L1SingletonCoeffs_subset, 
                             weights = 1/robust.se,
                             family = quasibinomial)
summary(LM_All_Interact_binom)
mean(1/L1SingletonCoeffs_subset$robust.se)

# Smooth proportions of 
PropSmoothed_InsL <- supsmu(L1SingletonCoeffs_subset$InsLength,
                            1*L1SingletonCoeffs_subset$blnSelect,
                            wt = 1/L1SingletonCoeffs_subset$robust.se)
plot(PropSmoothed_InsL$x, PropSmoothed_InsL$y, xlab = "Insertion length [bp]",
     ylab = "Proportion of L1 with selection signal", type = "l",
     ylim = c(0, 0.2), col = "red")
lines(PropSmoothed_InsL$x[PropSmoothed_InsL$x>=6000], 
      PropSmoothed_InsL$y[PropSmoothed_InsL$x>=6000])
InsLorder <- order(L1SingletonCoeffs_subset$InsLength)
lines(L1SingletonCoeffs_subset$InsLength[InsLorder], 
      LM_All_Interact_binom$fitted.values[InsLorder], col = "green")
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
                             weights = 1/L1SingletonCoeffs$robust.se, family  = quasibinomial)
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
L1SingletonCoeffs$FullDist2Loop <- L1SingletonCoeffs$Dist2Loop * L1SingletonCoeffs$blnFull
L1SingletonCoeffs$Dist2Domain <- Dist2Closest(L1SingletonCoeffs_GR, 
                                              DomainGR_Intersect)
L1SingletonCoeffs$Dist2LoopBins <- cut(L1SingletonCoeffs$Dist2Loop,
                                       breaks = seq(0, 24*10^6, 10^5))
max(L1SingletonCoeffs$Dist2Loop)

#######
# Regress against distance to loops
#######

# Polynomial regression coefficient against distance to loops for full L1
LM_Full1Mpp <- lm(coef ~ Dist2Loop + blnFull*Dist2Loop, data = L1SingletonCoeffs, 
               subset = robust.se > 0 & Dist2Loop <= 10^6,
               weights = 1/L1SingletonCoeffs$robust.se)
summary(LM_Full1Mpp)
LM_Full1 <- lm(coef ~ Dist2Loop, data = L1SingletonCoeffs, 
               subset = L1SingletonCoeffs$InsLength >= 6000 & L1SingletonCoeffs$robust.se > 0,
               weights = 1/L1SingletonCoeffs$robust.se)
LM_Full2 <- lm(coef ~ poly(Dist2Loop, 2), data = L1SingletonCoeffs, 
               subset = L1SingletonCoeffs$InsLength >= 6000 & L1SingletonCoeffs$robust.se > 0,
               weights = 1/L1SingletonCoeffs$robust.se)
LM_Full3 <- lm(coef ~ poly(Dist2Loop, 3), data = L1SingletonCoeffs, 
              subset = L1SingletonCoeffs$InsLength >= 6000 & L1SingletonCoeffs$robust.se > 0,
              weights = 1/L1SingletonCoeffs$robust.se)
AIC(LM_Full1)
AIC(LM_Full2)
AIC(LM_Full3)
LM_Full3$model
xVals <- L1SingletonCoeffs$Dist2Loop[L1SingletonCoeffs$InsLength >= 6000 & 
                                  L1SingletonCoeffs$robust.se > 0]
xValOrder <- order(xVals)
YVals <- LM_Full3$fitted.values[xValOrder]
xValsOrder <- xVals[xValOrder]
blnFull <- L1SingletonCoeffs$InsLength >= 6000
plot(L1SingletonCoeffs$Dist2Loop[blnFull], L1SingletonCoeffs$coef[blnFull],
     xlim = c(0, 10^6))
lines(xValsOrder, YVals)

# Polynomial regression coefficient against distance to loops for fragment L1
LM_Fragm1 <- lm(coef ~ Dist2Loop, data = L1SingletonCoeffs, 
               subset = L1SingletonCoeffs$InsLength < 6000 & L1SingletonCoeffs$robust.se > 0,
               weights = 1/L1SingletonCoeffs$robust.se)
LM_Fragm2 <- lm(coef ~ poly(Dist2Loop, 2), data = L1SingletonCoeffs, 
               subset = L1SingletonCoeffs$InsLength < 6000 & L1SingletonCoeffs$robust.se > 0,
               weights = 1/L1SingletonCoeffs$robust.se)
LM_Fragm3 <- lm(coef ~ poly(Dist2Loop, 3), data = L1SingletonCoeffs, 
               subset = L1SingletonCoeffs$InsLength < 6000 & L1SingletonCoeffs$robust.se > 0,
               weights = 1/L1SingletonCoeffs$robust.se)
AIC(LM_Fragm1)
AIC(LM_Fragm2)
AIC(LM_Fragm3)


LM_Fragm <- lm(coef ~ poly(Dist2Loop, 3), data = L1SingletonCoeffs, 
               subset = L1SingletonCoeffs$InsLength < 5900 & L1SingletonCoeffs$robust.se > 0,
               weights = 1/L1SingletonCoeffs$robust.se)
summary(LM_Fragm)
LM_Fragm$
LM_All <- lm(coef ~ Dist2Loop, data = L1SingletonCoeffs, 
               subset = L1SingletonCoeffs$robust.se > 0,
               weights = 1/L1SingletonCoeffs$robust.se)
summary(LM_All)

LM_All_Interact <- lm(coef ~ 
                        blnFull*poly(Dist2Loop, 3), data = L1SingletonCoeffs, 
             subset = L1SingletonCoeffs$robust.se > 0,
             weights = 1/L1SingletonCoeffs$robust.se)
summary(LM_All_Interact)

# Plot smoothed coefficients
CoefSmoothed_Full <- supsmu(L1SingletonCoeffs$Dist2Loop[L1SingletonCoeffs$blnFull],
                            1*L1SingletonCoeffs$coef[L1SingletonCoeffs$blnFull],
                            wt = 1/L1SingletonCoeffs$robust.se[L1SingletonCoeffs$blnFull])
CoefSmoothed_Fragm <- supsmu(L1SingletonCoeffs$Dist2Loop[!L1SingletonCoeffs$blnFull],
                             1*L1SingletonCoeffs$coef[!L1SingletonCoeffs$blnFull],
                             wt = 1/L1SingletonCoeffs$robust.se[!L1SingletonCoeffs$blnFull])
plot(L1SingletonCoeffs$Dist2Loop[L1SingletonCoeffs$blnFull],
     L1SingletonCoeffs$coef[L1SingletonCoeffs$blnFull], 
     xlab = "Distance to closest loop [bp]",
     ylab = "Coefficient", ylim = c(-1, 2))
points(L1SingletonCoeffs$Dist2Loop[L1SingletonCoeffs$blnFull],
     L1SingletonCoeffs$coef[L1SingletonCoeffs$blnFull], col = "red")
lines(CoefSmoothed_Fragm$x, CoefSmoothed_Fragm$y)
lines(CoefSmoothed_Full$x, CoefSmoothed_Full$y, col = "red")



Wts <- mean(L1SingletonCoeffs$robust.se[!is.na(blnSelect)])/
  L1SingletonCoeffs$robust.se
Wts[Wts == Inf] <- NA
mean(Wts, na.rm = T)

LM_All_Interact_binom <- glm(1*blnSelect ~ blnFull + blnFull*Dist2Loop, 
                             data = L1SingletonCoeffs_subset, 
                      subset = robust.se > 0 & (!is.na(blnSelect)),
                      weights = 1/robust.se, 
                      family  = quasibinomial)
LM_All_Interact_binom <- glm(1*blnSelect ~ 
                               blnFull*poly(Dist2Loop, 3), 
                             data = L1SingletonCoeffs_subset, 
                             subset = robust.se > 0 & (!is.na(blnSelect)),
                             weights = 1/robust.se/mean(1/robust.se), 
                             family  = quasibinomial)
summary(LM_All_Interact_binom)

mean(L1SingletonCoeffs)

LM_Full_binom <- glm(1*blnSelect ~ poly(Dist2Loop, 3), 
                    data = L1SingletonCoeffs, 
                    subset = robust.se > 0 & (!is.na(blnSelect)) & blnFull,
                    weights = 1/L1SingletonCoeffs$robust.se, 
                    family  = binomial)
sum(L1SingletonCoeffs$blnFull, na.rm = T)
summary(LM_All_binom)

LM_All_binom <- glm(1*blnSelect ~ Dist2Loop, 
                             data = L1SingletonCoeffs, 
                             subset = robust.se > 0 & (!is.na(blnSelect)),
                             weights = 1/L1SingletonCoeffs$robust.se, 
                             family  = quasibinomial)
summary(LM_All_binom)


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


# Polynomial regression coefficient against distance to loops for full L1
LM_Full1 <- glm(blnSelect ~ Dist2Loop, data = L1SingletonCoeffs, 
               subset = L1SingletonCoeffs$InsLength >= 6000 & L1SingletonCoeffs$robust.se > 0,
               weights = 1/L1SingletonCoeffs$robust.se,
               family = quasibinomial)
LM_Full2 <- glm(blnSelect ~ poly(Dist2Loop, 2), data = L1SingletonCoeffs, 
               subset = L1SingletonCoeffs$InsLength >= 6000 & L1SingletonCoeffs$robust.se > 0,
               weights = 1/L1SingletonCoeffs$robust.se,
               family = quasibinomial)
LM_Full3 <- glm(blnSelect ~ poly(Dist2Loop, 3), data = L1SingletonCoeffs, 
               subset = L1SingletonCoeffs$InsLength >= 6000 & L1SingletonCoeffs$robust.se > 0,
               weights = 1/L1SingletonCoeffs$robust.se,
               family = quasibinomial)
LM_Full1$
LM_Full2$aic
LM_Full3$aic
summary(LM_Full2)


# Polynomial regression coefficient against distance to loops for fragment L1
LM_Fragm1 <- glm(blnSelect ~ Dist2Loop, data = L1SingletonCoeffs, 
                subset = L1SingletonCoeffs$InsLength < 6000 & L1SingletonCoeffs$robust.se > 0,
                weights = 1/L1SingletonCoeffs$robust.se,
                family = quasibinomial)
LM_Fragm2 <- glm(blnSelect ~ poly(Dist2Loop, 2), data = L1SingletonCoeffs, 
                subset = L1SingletonCoeffs$InsLength < 6000 & L1SingletonCoeffs$robust.se > 0,
                weights = 1/L1SingletonCoeffs$robust.se,
                family = quasibinomial)
LM_Fragm3 <- glm(blnSelect ~ poly(Dist2Loop, 3), data = L1SingletonCoeffs, 
                subset = L1SingletonCoeffs$InsLength < 6000 & L1SingletonCoeffs$robust.se > 0,
                weights = 1/L1SingletonCoeffs$robust.se,
                family = quasibinomial)
summary(LM_Fragm1)
summary(LM_Fragm2)
summary(LM_Fragm3)


##########################
#                        #
#    Compare distances   #
#                        #
##########################

Distances  <- L1SingletonCoeffs$Dist2Loop
blnSelect  <- L1SingletonCoeffs$blnSelect
blnFull    <- L1SingletonCoeffs$blnFull
SampleSize = 10000
#AnalyzeDistanceInteraction <- function(Distances, blnSelect, 
#   blnFull, SampleSize = 10000){
blnNA           <- is.na(Distances) | is.na(blnSelect) | is.na(blnFull)
Distances       <- Distances[!blnNA]
blnSelect       <- blnSelect[!blnNA]
blnFull         <- blnFull[!blnNA]

blnSelectFull   <- blnSelect[blnFull]
blnSelectFragm  <- blnSelect[!blnFull]
Distances_Full  <- Distances[blnFull]
Distances_Fragm <- Distances[!blnFull]

MeanDiffFull <- mean(Distances_Full[blnSelectFull]) - 
  mean(Distances_Full[!blnSelectFull])
MeanDiffFragm <- mean(Distances_Fragm[blnSelectFragm]) - 
  mean(Distances_Fragm[!blnSelectFragm])
MeanDiffBoth  <- MeanDiffFull - MeanDiffFragm
NrFull        <- sum(blnFull)
NrFragm       <- sum(!blnFull)
NrSelectFull  <- sum(blnSelectFull)
NrSelectFragm <- sum(blnSelectFragm)


SampledDiffMat <- sapply(1:SampleSize, function(x) {
  SampleSelectFull   <- sample(1:NrFull,  size = NrSelect)
  SampleSelectFragm  <- sample(1:NrFragm, size = NrSelectFragm)
  MeanDiffFull <- mean(Distances_Full[SampleSelectFull]) - 
    mean(Distances_Full[-SampleSelectFull])
  MeanDiffFragm <- mean(Distances_Fragm[SampleSelectFragm]) - 
    mean(Distances_Fragm[-SampleSelectFragm])
  MeanDiffBoth = MeanDiffFull - MeanDiffFragm
  c(MeanDiffFull = MeanDiffFull, MeanDiffFragm = MeanDiffFragm,
    MeanDiffBoth = MeanDiffBoth)
  
})

2*sum(SampledDiffMat["MeanDiffFull", ]  <= MeanDiffFull) / SampleSize
2*sum(SampledDiffMat["MeanDiffFragm", ] >= MeanDiffFragm) / SampleSize
2*sum(SampledDiffMat["MeanDiffBoth", ]  <= MeanDiffBoth) / SampleSize

fisher.test(rbind(c(NrSelectFull, NrFull - NrSelectFull),
                  c(NrSelectFragm, NrFragm - NrSelectFragm)))

##########################
#                        #
#    Save image          #
#                        #
##########################


# save.image(OutputPath)

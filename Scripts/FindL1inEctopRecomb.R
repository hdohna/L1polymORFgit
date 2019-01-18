# The script below identifies deletions that are likely L1 polymorphis,s

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(GenomicRanges)
library(rtracklayer)

# Specify file paths
DataFolder          <- 'D:/L1polymORF/Data/'
L1GRPath            <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
L1RefRangePath      <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
L1RefRangeBedPath   <- 'D:/L1polymORF/Data/L1allRefRanges_hg19.bed'

# Read in ranges of non-polymorphic
# L1NoPolyGR <- import.bed("D:/L1polymORF/Data/L1GRangesNoPoly.bed")
# load(L1GRPath)
# load(L1RefRangePath)
# L1GRangesAll <- c(L1GRanges, L1_1000G_GR_hg19)
# L1Width <- c(width(L1GRanges), 
#              as.numeric(as.character(L1_1000G$L1End)) - 
#                as.numeric(as.character(L1_1000G$L1Start)))
# L1GRangesAll <- L1GRangesAll[!duplicated(start(L1GRangesAll))]
# L1Width <- L1Width[!duplicated(start(L1GRangesAll))]

# Read in ranges of all L1s (L1HS, L1PA, etc.)
L1GRangesAll <- import.bed('D:/L1polymORF/Data/L1allRefRanges_hg19.bed')

# Get L1 table
L1Table <- read.csv('D:/L1polymORF/Data/L1Table.csv', as.is = T)
L1Table <- L1Table[nchar(L1Table$chromosome) <= 2, ]
PosCharTable <- paste(L1Table$chromosome, L1Table$genoStart, L1Table$genoEnd)
PosCharGR <- paste(as.character(seqnames(L1GRangesAll)), 
                   start(L1GRangesAll), end(L1GRangesAll))
PosMatch <- match(PosCharGR, PosCharTable)
L1Table <- L1Table[PosMatch, ]

# Subset GRanges to get only L1PA 
blnL1PA <- substr(L1Table$repName, 1, 4) == "L1PA"
L1GRanges_L1PA <- L1GRangesAll[blnL1PA]

# Get file names, loop over files and do the filtering
AllFiles <- list.files(DataFolder, pattern = "_L1NoPoly_DelDupInv", 
                       full.names = T)
AllFiles <- list.files(DataFolder, pattern = "_L1all_DelDupInv", 
                       full.names = T)
StartV <- NULL
EndV   <- NULL
ChromV <- NULL
InfoV  <- NULL
InfoV_all  <- NULL
LineCount <- 0
for (InFile in AllFiles){
  VcfFile <- read.table(InFile, stringsAsFactors = F)
  LineCount <- LineCount + nrow(VcfFile)
  Chrom   <- strsplit(InFile, "_")[[1]][2]
  cat("Analyzing", Chrom, ":\n")
  cat("Total number of entries", LineCount, ":\n")
  idxCNV  <- which(VcfFile$V5 == "<CN0>")
  StartV  <- c(StartV, VcfFile$V2[idxCNV])
  ChromV  <- c(ChromV, rep(Chrom, length(idxCNV)))
  InfoV_all   <- c(InfoV, VcfFile$V8)
  InfoV   <- c(InfoV, VcfFile$V8[idxCNV])
  Ends <- sapply(VcfFile$V8[idxCNV], function(x){
    Split1 <- strsplit(x, ";")[[1]]
    EndPart <- grep("END=", Split1, value = T)
    if(length(grep("CIEND=", EndPart)) > 0){
      EndPart <- EndPart[-grep("CIEND=", EndPart)]
    }
    if (length(EndPart) > 0){
      as.numeric(strsplit(EndPart, "END=")[[1]][2])
    } else {
      NA
    }
    
  })
  EndV <- c(EndV, unlist(Ends))
  cat("Added", length(idxCNV), "start and chromosome values and", length(Ends),
      "end values\n\n")
  
}

length(EndV)

# Create genomic ranges for start and end of variants
StartGR <- GRanges(seqnames = substr(ChromV, 4, nchar(ChromV)),
                   IRanges(start = StartV, end = StartV))
EndGR   <- GRanges(seqnames = substr(ChromV, 4, nchar(ChromV)),
                   IRanges(start = EndV, end = EndV))
GrChar <- paste(ChromV, StartV, EndV, sep = "_")
StartGR <- StartGR[!duplicated(GrChar)]
EndGR   <- EndGR[!duplicated(GrChar)]

# StartGR <- GRanges(seqnames = ChromV, 
#                    IRanges(start = StartV - 10, end = StartV + 10))
# EndGR   <- GRanges(seqnames = ChromV, 
#                    IRanges(start = EndV - 10, end = EndV + 10))

# Remove SVs that overlap with the same L1 on start and end
OL_Start <- findOverlaps(StartGR, L1GRangesAll)
OL_End   <- findOverlaps(EndGR, L1GRangesAll)
matchStartEnd <- match(OL_Start@from, OL_End@from)
idxSameL1 <- which(OL_Start@to == OL_End@to[matchStartEnd])
StartGR <- StartGR[-idxSameL1]
EndGR   <- EndGR[-idxSameL1]

# Get the L1 5' position of the start and end of the start and end of
# L1
OL_Start <- findOverlaps(StartGR, L1GRangesAll)
OL_End   <- findOverlaps(EndGR, L1GRangesAll)
StartL15P <- rep(NA, length(StartGR))
StartL15P_n <- sapply(1:length(OL_Start@to), function(x){
  if (L1Table$strand[OL_Start@to[x]] == "+"){
    start(StartGR)[OL_Start@from[x]] - start(L1GRangesAll)[OL_Start@to[x]]
  } else {
    end(L1GRangesAll)[OL_Start@to[x]] - start(StartGR)[OL_Start@from[x]]
  }
})
StartL15P[OL_Start@from] <- StartL15P_n
EndL15P <- rep(NA, length(EndGR))
EndL15P_n <- sapply(1:length(OL_End@to), function(x){
  if (L1Table$strand[OL_End@to[x]] == "+"){
    start(EndGR)[OL_End@from[x]] - start(L1GRangesAll)[OL_End@to[x]]
  } else {
    end(L1GRangesAll)[OL_End@to[x]] - start(EndGR)[OL_End@from[x]]
  }
})
EndL15P[OL_End@from] <- EndL15P_n
DiffV <- seq(0, 200, 10)
NrSamePos <- sapply(DiffV, function(x) {
  sum(abs(StartL15P - EndL15P) <= x, na.rm = T)
})
plot(DiffV, NrSamePos)
hist(abs(StartL15P - EndL15P))

# Get the L1 name
StartL1Name <- rep(NA, length(StartGR))
table(L1Table$repName)
StartL1Name[OL_Start@from] <- L1Table$repName[OL_Start@to]
EndL1Name <- rep(NA, length(EndGR))
EndL1Name[OL_End@from] <- L1Table$repName[OL_End@to]
idxSameName <- which(StartL1Name == EndL1Name)
sum(StartL1Name == EndL1Name, na.rm = T)
sample.int(5)

# Randomly reorder start name and determine the number of 
idxNotNA <- which((!is.na(StartL1Name)) & (!is.na(EndL1Name)))
StartL1Name_NoNA <- StartL1Name[idxNotNA]
EndL1Name_NoNA   <- EndL1Name[idxNotNA]
L1Combo <- paste(StartL1Name_NoNA, EndL1Name_NoNA)
L1ComboCount <- table(L1Combo)
sort(L1ComboCount, decreasing = T)[1:50]
L1ComboCount["L1HS L1HS"]
sum(StartL1Name_NoNA == EndL1Name_NoNA)
SampledL1NameMatch <- sapply(1:10000, function(z){
  idx <- sample.int(length(EndL1Name_NoNA))
  sum(StartL1Name_NoNA == EndL1Name_NoNA[idx])
})
hist(SampledL1NameMatch, xlim = c(0, 800))

# Create a predictor variable for involvement in ectopic recombination
L1Width        <- width(L1GRangesAll)
L1Width[L1Width >= 6200] <- 6200
L1WidthOrder   <- order(L1Width, decreasing = T)
OrderMatch     <- match(1:length(L1GRangesAll), L1WidthOrder)
Deltas         <- c(L1Width[L1WidthOrder[1]], 
                    L1Width[L1WidthOrder[-length(L1Width)]] - L1Width[L1WidthOrder[-1]])
DeltasSqProd   <- 10^-10*Deltas^2 * (L1WidthOrder - 1)
Rev            <- length(L1GRangesAll):1
L1WidthProd    <- cumsum(DeltasSqProd[Rev])[Rev]
RecPredict     <- L1WidthProd[OrderMatch]

# Create a predictor variable for involvement in ectopic recombination
# for L1PA only
L1Width_L1PA        <- width(L1GRanges_L1PA)
L1Width_L1PA[L1Width_L1PA >= 4500] <- 4500
L1Width_L1PAOrder   <- order(L1Width_L1PA, decreasing = T)
OrderMatch     <- match(1:length(L1GRanges_L1PA), L1Width_L1PAOrder)
Deltas         <- c(L1Width_L1PA[L1Width_L1PAOrder[1]], 
                    L1Width_L1PA[L1Width_L1PAOrder[-length(L1Width_L1PA)]] - 
                      L1Width_L1PA[L1Width_L1PAOrder[-1]])
DeltasSqProd   <- 10^-10*Deltas^2 * (L1Width_L1PAOrder - 1)
Rev            <- length(L1GRanges_L1PA):1
L1Width_L1PAProd    <- cumsum(DeltasSqProd[Rev])[Rev]
RecPredict_L1PA     <- L1Width_L1PAProd[OrderMatch]

# Get indicator for L1s that overlap with either end of the CNV
# blnOL   <- overlapsAny(L1NoPolyGR, StartGR) | overlapsAny(L1NoPolyGR, EndGR)
blnOL_Start <- overlapsAny(StartGR, L1GRangesAll) 
blnOL_End   <- overlapsAny(EndGR, L1GRangesAll) 
blnOL_Both  <- blnOL_Start & blnOL_End
blnOL       <- overlapsAny(L1GRangesAll, StartGR) | overlapsAny(L1GRangesAll, EndGR)
OLCount     <- countOverlaps(L1GRangesAll, StartGR) + countOverlaps(L1GRangesAll, EndGR)
OLCount2    <- countOverlaps(L1GRangesAll, StartGR[blnOL_Both]) + 
               countOverlaps(L1GRangesAll, EndGR[blnOL_Both])
OLCount3    <- countOverlaps(L1GRangesAll, StartGR[idxSameName]) + 
  countOverlaps(L1GRangesAll, EndGR[idxSameName])

blnOLBoth    <- overlapsAny(L1GRangesAll, StartGR[blnOL_Both]) | 
  overlapsAny(L1GRangesAll, EndGR[blnOL_Both])
PropOLStart <- sum(blnOL_Start) / length(StartGR)
PropOLEnd   <- sum(blnOL_End) / length(StartGR)
PropOLBoth  <- sum(blnOL_Both) / length(StartGR)
pbinom(sum(blnOL_Both), length(StartGR), PropOLStart*PropOLEnd)

# Read in file and create GRanges
# RecData <- read.delim("D:/L1polymORF/Data/hg19deCodeRecomb.txt")
RecData <- read.delim("D:/L1polymORF/Data/hg19RecombRate.txt", 
                      stringsAsFactors = F)
RecData$chrom <- substr(RecData$chrom, 4, nchar(RecData$chrom))
Rec_GR <- makeGRangesFromDataFrame(RecData)

# Find overlaps to catalog L1 and create a vector of recombination values
cat("Calculating recombination rate per L1")
RecL1CatOverlaps <- findOverlaps(L1GRangesAll, Rec_GR)
MeanRec         <- aggregate(RecData$decodeAvg[RecL1CatOverlaps@to], 
                             by = list(RecL1CatOverlaps@from) , FUN = mean) 
MeanRecPerL1 <- rep(NA, length(L1GRangesAll))
MeanRecPerL1[MeanRec$Group.1] <- MeanRec$x

# Create a vector of L1 width classes
L1widthClass <- cut(width(L1GRangesAll), breaks = 
                      seq(0, 6250, 250))
L1widthClass_L1PA <- cut(width(L1GRanges_L1PA), breaks = 
                      seq(0, 6250, 250))

# Get mean L1 overlap with deletion per width class
L1widthAggregated <- aggregate(data.frame(L1Width = width(L1GRangesAll), 
                                          blnOL = blnOL,
                                          blnOLBoth = blnOLBoth,
                                          OLCount = OLCount,
                                          OLCount2 = OLCount2,
                                          OLCount3 = OLCount3), 
                               by = list(L1widthClass), 
                               FUN = mean)
L1widthAggregated_var <- aggregate(data.frame(L1Width = width(L1GRangesAll), 
                                          blnOL = blnOL,
                                          blnOLBoth = blnOLBoth,
                                          OLCount = OLCount,
                                          OLCount2 = OLCount2,
                                          OLCount3 = OLCount3), 
                               by = list(L1widthClass), 
                               FUN = var)
L1widthAggregated_n <- aggregate(data.frame(L1Width = width(L1GRangesAll), 
                                              blnOL = blnOL,
                                              blnOLBoth = blnOLBoth,
                                              OLCount = OLCount,
                                              OLCount2 = OLCount2,
                                              OLCount3 = OLCount3), 
                                   by = list(L1widthClass), 
                                   FUN = length)
L1widthAgg_L1PA <- aggregate(data.frame(L1Width = width(L1GRanges_L1PA), 
                                          blnOL = blnOL[blnL1PA],
                                          blnOLBoth = blnOLBoth[blnL1PA],
                                          OLCount = OLCount[blnL1PA],
                                          OLCount2 = OLCount2[blnL1PA],
                                          OLCount3 = OLCount3[blnL1PA]), 
                               by = list(L1widthClass_L1PA), 
                               FUN = mean)
L1widthAgg_L1PA_var <- aggregate(data.frame(L1Width = width(L1GRanges_L1PA), 
                                            blnOL = blnOL[blnL1PA],
                                            blnOLBoth = blnOLBoth[blnL1PA],
                                            OLCount = OLCount[blnL1PA],
                                            OLCount2 = OLCount2[blnL1PA],
                                            OLCount3 = OLCount3[blnL1PA]), 
                                 by = list(L1widthClass_L1PA), 
                                   FUN = var)
L1widthAgg_L1PA_n <- aggregate(data.frame(L1Width = width(L1GRanges_L1PA), 
                                            blnOL = blnOL[blnL1PA],
                                            blnOLBoth = blnOLBoth[blnL1PA],
                                            OLCount = OLCount[blnL1PA],
                                            OLCount2 = OLCount2[blnL1PA],
                                            OLCount3 = OLCount3[blnL1PA]), 
                                 by = list(L1widthClass_L1PA), 
                                 FUN = length)

# Check whether probability of overlap depends on end and full-length indicator
blnFull <- width(L1GRangesAll)  >= 4700
GLM_OL <- glm(OLCount ~ width(L1GRangesAll) + blnFull, family = poisson("identity"))
summary(GLM_OL)
GLM_OL_log <- glm(OLCount ~ log(width(L1GRangesAll)) + blnFull, family = poisson)
summary(GLM_OL_log)

GLM_OL2_lin <- glm(OLCount2 ~ width(L1GRangesAll) + blnFull, family = poisson("identity"),
               start = c(0, 1.2*10^-5, 0))
summary(GLM_OL2_lin)

GLM_OL2_log <- glm(OLCount2 ~ log(width(L1GRangesAll)) + blnFull, family = poisson,
                   start = c(0, 1.2*10^-5, 0))
summary(GLM_OL2_log)

GLM_OL3_lin <- glm(OLCount3 ~ width(L1GRangesAll) + blnFull, 
                   family = poisson("identity"),
               start = c(0, 7*10^-6, 0))
summary(GLM_OL3_lin)
GLM_OL3_log <- glm(OLCount3 ~ log(width(L1GRangesAll)) + blnFull, family = poisson,
               start = c(0, 7*10^-6, 0))
summary(GLM_OL3_log)

GLM_OL3_pred <- glm(OLCount3 ~ log(RecPredict + 10^-10) + blnFull, 
                    family = poisson,
                   start = c(1.4, 1, -5))
summary(GLM_OL3_pred)
# GLM_OL3_pred_ident <- glm(OLCount3 ~ RecPredict +  MeanRecPerL1,
#                           family = poisson("identity"),
#                           start = c(10^-5, 1.4, 0))
GLM_OL3_pred_ident <- glm(OLCount3 ~ 0 + RecPredict + blnFull,
                          family = poisson("identity"),
                          start = c(10^-5, 1.4, 0))
summary(GLM_OL3_pred_ident)
GLM_OL3_pred_ident_L1PA <- glm(OLCount3[blnL1PA] ~ RecPredict_L1PA,
                    family = poisson("identity"),
                    start = c(10^-5, 1.4))
blnNARec <- is.na(MeanRecPerL1)

# Get order of L1 width
WidthOrder <- order(width(L1GRangesAll))
WidthOrderRec <- order(width(L1GRangesAll)[!blnNARec])
WidthOrder_L1PA <- order(width(L1GRanges_L1PA))

# Get a smoothed curve through deletions as function of Length
DelVsL1Length <- supsmu(x = width(L1GRangesAll), y = OLCount3)
save(list = "DelVsL1Length", file = "D:/L1polymORF/Data/DelVsL1Length.RData")
xOrder <- order(DelVsL1Length$x)
plot(L1widthAggregated$L1Width, L1widthAggregated$OLCount, xlab = "L1 length",
     ylab = "Proportion of deletions with end in L1", ylim = c(0, 0.14))
lines(width(L1GRangesAll)[WidthOrder], GLM_OL$fitted.values[WidthOrder])

plot(L1widthAggregated$L1Width, L1widthAggregated$OLCount2, 
     xlab = "L1 length", ylab = "Proportion of deletions flanked by L1")
lines(width(L1GRangesAll)[WidthOrder], GLM_OL2_lin$fitted.values[WidthOrder],
      col = "red")
lines(width(L1GRangesAll)[WidthOrder], GLM_OL2_log$fitted.values[WidthOrder],
      col = "blue")

plot(L1widthAggregated$L1Width, L1widthAggregated$OLCount3, 
     xlab = "L1 length", ylab = "Deletions per L1",
     ylim = c(0, 0.06))
AddErrorBars(MidX = L1widthAggregated$L1Width, 
             MidY = L1widthAggregated$OLCount3, 
             ErrorRange = sqrt(L1widthAggregated_var$OLCount3 /
                                 L1widthAggregated_n$OLCount3),
             TipWidth = 20)
lines(DelVsL1Length$x[xOrder], DelVsL1Length$y[xOrder])
lines(width(L1GRangesAll)[WidthOrder], GLM_OL3_pred$fitted.values[WidthOrder],
      col = "red")
lines(width(L1GRangesAll)[WidthOrder], 
      GLM_OL3_pred_ident$fitted.values[WidthOrder], col = "red")
# lines(width(L1GRangesAll)[!blnNARec][WidthOrderRec], 
#       GLM_OL3_pred_ident$fitted.values[WidthOrderRec], col = "red")
length(GLM_OL3_pred_ident$fitted.values)
length(L1GRangesAll)
# lines(width(L1GRangesAll)[WidthOrder], 1.4*RecPredict[WidthOrder],
#       col = "red")
# lines(width(L1GRangesAll)[WidthOrder], GLM_OL3_lin$fitted.values[WidthOrder],
#       col = "red")
lines(width(L1GRangesAll)[WidthOrder], GLM_OL3_log$fitted.values[WidthOrder])

# Plot mean deletions per L1 size class for L1PA
plot(L1widthAgg_L1PA$L1Width, L1widthAgg_L1PA$OLCount3, 
     xlab = "L1 length", ylab = "Deletions per L1",
     ylim = c(0, 0.06))
AddErrorBars(MidX = L1widthAgg_L1PA$L1Width, 
             MidY = L1widthAgg_L1PA$OLCount3, 
             ErrorRange = sqrt(L1widthAgg_L1PA_var$OLCount3 /
                                 L1widthAgg_L1PA_n$OLCount3),
             TipWidth = 20)
lines(width(L1GRanges_L1PA)[WidthOrder_L1PA], 
      GLM_OL3_pred_ident_L1PA$fitted.values[WidthOrder_L1PA],
      col = "red")

# Sample regression coefficients
NrSamples <- 1000
# SampledCoeffs <- sapply(1:NrSamples, function(x){
#   SampledL1 <- sample.int(length(L1GRangesAll), size = sum(blnOL),
#                           prob = width(L1GRangesAll))
#   blnOL_local <- c(1:length(L1GRangesAll)) %in% SampledL1
#   GLM_OL2 <- glm(blnOL_local ~ width(L1GRangesAll) + blnFull, family = binomial)
#   summary(GLM_OL2)$coefficients[2:3,'Estimate']
#   
# })
# dim(SampledCoeffs)
# sum(SampledCoeffs[1, ] >= summary(GLM_OL2)$coefficients[2,'Estimate'] ) / NrSamples
# sum(SampledCoeffs[2, ] >= summary(GLM_OL2)$coefficients[3,'Estimate'] ) / NrSamples
# hist(SampledCoeffs[1, ])
# hist(SampledCoeffs[2, ])

# Get type of SV
SVtype <- sapply(InfoV, function(x){
  Split1 <- strsplit(x, ";")[[1]]
  grep("SVTYPE=", Split1, value = T)
})
table(SVtype)

# Read in SVcounts
SVcounts <- read.table("D:/L1polymORF/Data/SVcountAllChroms")
sum(blnOL_Both) / sum(SVcounts$V1)

# Read chromosome lengths 
load('D:/L1polymORF/Data/ChromLengthsHg19.Rdata')
PropL1 <- sum(width(L1GRangesAll)) / sum(ChromLengthsHg19)
PropL1^2
pbinom(sum(blnOL_Both), length(StartGR), PropL1^2, lower.tail = F)
pbinom(sum(blnOL_Start), length(StartGR), PropL1, lower.tail = F)
pbinom(sum(blnOL_End), length(StartGR), PropL1, lower.tail = F)
PropOLStart
PropOLEnd

# Test Poisson regression
xVals <- seq(0.1, 5, 0.1)
yVals <- rpois(n = length(xVals), lambda = xVals)
plot(xVals, yVals)

GL1 <- glm(yVals ~ xVals, family = poisson)
summary(GL1)
GL2 <- glm(yVals ~ xVals, family = poisson("identity"),
           start = c(0, 1))
summary(GL2)
GL3 <- glm(yVals ~ log(xVals), family = poisson,
           start = c(0, 1))
summary(GL3)

# The following script analyzes the output creates by isPcr

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.R')

# load packages
library(seqinr)
library(ape)
library(rtracklayer)

#########################################
#                                       #
#           Read in data                #
#                                       #
#########################################

# Specify file path to isPcr output
isPcrOutputFolder <- "D:/L1polymORF/Data/L1Primers"
isPcrOutputPrefix <- "L1catalog.isPcr.output."

# Read primer 3 info
Info_P3 <- read.csv("D:/L1polymORF/Data/L1catalog.primer3.info.csv", as.is = T,
                    stringsAsFactors = F, row.names = NULL)
Info_P3$isLeft <- sapply(Info_P3$Name, function(x){
 strsplit(x, "_")[[1]][2] == "S1"
})
  
# Read output and create an info table
isPcrOutputFiles <- list.files(isPcrOutputFolder, pattern = isPcrOutputPrefix,
                               full.names = T)
Info_isPCR <- data.frame()
for (isPcrOutputPath in isPcrOutputFiles){
  
  isPcrOutput <- read.dna(isPcrOutputPath, format = "fasta", as.character = T)
  Info_isPCRPerChr <- read.table(text = names(isPcrOutput), 
                     col.names = c("Position", "Name", "Length","Primer1", "Primer2"),
                     as.is = T)
  Info_isPCR <- rbind(Info_isPCR, Info_isPCRPerChr)
}

# Append chromosome and accession number
Info_isPCR$Chromosome <- sapply(Info_isPCR$Position, function(x) strsplit(x, ":")[[1]][1])
Info_isPCR$Accession  <- sapply(Info_isPCR$Name, function(x) strsplit(x, "_")[[1]][1])
Info_isPCR$LengthNum  <- as.numeric(substr(Info_isPCR$Length, 1, nchar(Info_isPCR$Length) - 2))
Info_isPCR$start_HG38 <- sapply(Info_isPCR$Position, function(x) {
                          Split1 <- strsplit(x, ":")[[1]][2]
                          if (length(grep("\\+", Split1)) > 0){
                            as.numeric(strsplit(Split1, "\\+")[[1]][1])
                          } else {
                            as.numeric(strsplit(Split1, "\\-")[[1]][1])
                          }
                        })
Info_isPCR$end_HG38 <- sapply(Info_isPCR$Position, function(x) {
  Split1 <- strsplit(x, ":")[[1]][2]
  if (length(grep("\\+", Split1)) > 0){
    as.numeric(strsplit(Split1, "\\+")[[1]][2])
  } else {
    as.numeric(strsplit(Split1, "\\-")[[1]][2])
  }
})

# Read in L1 catalog 
L1Catalog <- read.csv("D:/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv",
                      as.is = T)

# Indicator whether L1 is in reference
blnInReference <- (L1Catalog$end_HG38 - L1Catalog$start_HG38) > 5000
AccNotInRef  <- L1Catalog$Accession[which(!blnInReference)]

#########################################
#                                       #
#       Define functions                #
#                                       #
#########################################

# Function to determine best primers for a particular L1
L1Name <- "AL022308_S1"
GetBestPrimer <- function(L1Name, Pr3Info = Info_P3){
  idxRows   <- grep(L1Name, Pr3Info$Name)
  Pr3subset <- Pr3Info[idxRows, ]
  MeanGC    <- rowMeans(Pr3subset[ ,c("GC_Left", "GC_Right")])
  OptDiffLeft  <- abs(50 - Pr3subset$GC_Left)
  OptDiffRight <- abs(50 - Pr3subset$GC_Right)
  OptDiff      <- 0.5*(OptDiffLeft + OptDiffRight) 
  Pr3subset <- Pr3subset[order(OptDiff), ]
  blnPass   <- 
    Pr3subset$GC_Left >= 40 & Pr3subset$GC_Left <= 60 &
    Pr3subset$GC_Right >= 40 & Pr3subset$GC_Right <= 60 &
    Pr3subset$Tm_Left >= 55 & Pr3subset$Tm_Left <= 65 &
    Pr3subset$Tm_Right >= 55 & Pr3subset$Tm_Right <= 65 &
    Pr3subset$ProductSize >= 400  &
    Pr3subset$Length_Left >= 20 & Pr3subset$Length_Right >= 20 
  blnAnyPassFilter <- any(blnPass)
  if (blnAnyPassFilter){
    OptRow <- min(which(blnPass))
  } else {
    warning("No primer passed filter, returning the one with best Tm")
    OptRow <- 1
  }
  data.frame(Pr3subset[OptRow, ], PassFilter = blnAnyPassFilter, stringsAsFactors = F)
}

#########################################
#                                       #
#           Subset Info_P3              #
#                                       #
#########################################

# Match accession numbers
Info_P3$Accession  <- sapply(Info_P3$Name, function(x) strsplit(x, "_")[[1]][1])
AccMatch <- match(Info_P3$Accession, L1Catalog$Accession)

# Loop over all primers in Info_P3 and retain only the ones whose alternative 
# product is greater than 600 bp or smaller than 300 bp
x <- Info_P3$Name[1]
i <- 1
InfoCompare <- t(sapply(1:nrow(Info_P3), function(i){
  x <- Info_P3$Name[i]
  Info_isPCRSubset <- Info_isPCR[Info_isPCR$Name == x, ]
  if (Info_P3$isLeft[i]){
    FocalRegion <- L1Catalog$start_HG38[AccMatch][i]
  } else {
    FocalRegion <- L1Catalog$end_HG38[AccMatch][i]
  }
  if (is.na(FocalRegion)){
    blnRegionCovered  <- L1Catalog$Chromosome[AccMatch][i] == Info_isPCRSubset$Chromosome
    blnInRef <- T
  } else {
    blnRegionCovered  <- FocalRegion > Info_isPCRSubset$start_HG38 &
      FocalRegion < Info_isPCRSubset$end_HG38
    blnInRef <- blnInReference[AccMatch][i]
  }
  SizeDiff    <- abs(Info_P3$ProductSize[i] - Info_isPCRSubset$LengthNum)
  blnSameSize <- SizeDiff < 100
  if (any(!blnRegionCovered)){
    idxNotCovered  <- which(!blnRegionCovered)
    idxCloseSize   <- which.min(SizeDiff[idxNotCovered])
    AltProductSize <- Info_isPCRSubset$LengthNum[idxNotCovered[idxCloseSize]]
  } else {
    AltProductSize <- NA
  }
  if (blnInRef){
    PrimerOK <- sum(blnRegionCovered & blnSameSize) == 1 &
      sum(!blnRegionCovered & blnSameSize) == 0
  } else {
    PrimerOK <- sum(blnSameSize) == 0
  }
  c(PrimerOK = PrimerOK,
    NrAlternativeProduct = sum(!blnRegionCovered),
    AltProductSize = AltProductSize)
}))
sum(InfoCompare[,"PrimerOK"])

# Subset Info_P3 to get the OK primer
blnP3subset   <- InfoCompare[,"PrimerOK"] == 1
Info_P3Subset <- Info_P3[blnP3subset, ]
Info_P3Subset$AltProductSize <- InfoCompare[blnP3subset,"AltProductSize"]

# Get all sequence names
SeqNames <- sapply(as.character(Info_P3Subset$Name), function(x){
  SplitName <- strsplit(x, "_")[[1]]
  paste(SplitName[1:2], collapse = "_")
})

# Select for each sequence name the best primer
BestPrimers <- data.frame(stringsAsFactors = F)
for (SeqName in unique(SeqNames)){
  BestPrimer  <- GetBestPrimer(SeqName, Pr3Info = Info_P3Subset)
  BestPrimers <- rbind(BestPrimers, BestPrimer)
}

# Determine for each accession in L1 how often it is covered in BestPrimers
AccCount <- sapply(unique(L1Catalog$Accession), function(x) sum(BestPrimers$Accession == x, na.rm = T))
table(AccCount)

# Add to each primer whether it is 5' or 3'
AccMatchBP <- match(BestPrimers$Accession, L1Catalog$Accession)
BestPrimers$is5Prime <- FALSE
BestPrimers$is5Prime <- BestPrimers$isLeft & L1Catalog$Strand[AccMatchBP] == "+"
BestPrimers$inRef    <- blnInReference[AccMatchBP]

# Select for each accession number the appropriate primers
AllAccs <- unique(BestPrimers$Accession)
BestPrimerSubset <- t(sapply(AllAccs, function(x){
  BPSubset    <- BestPrimers[BestPrimers$Accession == x, ]
  if (sum(BPSubset$PassFilter) < 0){
    BPSubset    <- BPSubset[BPSubset$PassFilter, ]
  }
  L1CatSubset <- L1Catalog[which(L1Catalog$Accession == x)[1], ]
  if (nrow(BPSubset) == 2){
    if (L1CatSubset$Strand == "+") {
      BPSubset <- BPSubset[BPSubset$isLeft, ]
    } else {
      BPSubset <- BPSubset[!BPSubset$isLeft, ]
    }
  }
  BPSubset
}))

BestPrimerSubset <- as.data.frame(BestPrimerSubset)

# Retain the primers that pass the filter
BestPrimerSubset <- BestPrimerSubset[unlist(BestPrimerSubset$PassFilter),]

# Retain primers for 96 L1 that are in reference
inRefVect <- unlist(BestPrimerSubset$inRef)
P5Vect    <- unlist(BestPrimerSubset$is5Prime)
NewOrder  <- order(inRefVect, !P5Vect)
BestPrimerSubset <- BestPrimerSubset[NewOrder[1:96], ]

# Reorder primer table
ColsLeft <- c("Name", "ProductSize", "PassFilter", 
              grep("_Left", colnames(BestPrimerSubset), value = T))
ColsRight <- c("Name", "ProductSize", "PassFilter", 
              grep("_Right", colnames(BestPrimerSubset), value = T))
BestPrimerLeft  <- as.data.frame(BestPrimerSubset[ ,ColsLeft])
BestPrimerRight <- as.data.frame(BestPrimerSubset[ ,ColsRight])
BestPrimerLeft$PrimerType  <- "Left"
BestPrimerRight$PrimerType <- "Right"
colnames(BestPrimerLeft)  <- gsub("_Left", "", colnames(BestPrimerLeft))
colnames(BestPrimerRight) <- gsub("_Right", "", colnames(BestPrimerRight))
BestPrimerReordered <- rbind(BestPrimerLeft, BestPrimerRight)
BestPrimerReordered <- as.data.frame(BestPrimerReordered)
NameOrder           <- order(as.character(BestPrimerReordered$Name))
BestPrimerReordered <- BestPrimerReordered[NameOrder, ]

# Determine how many non-reference 
L1Catalog$WithPrimers <- sapply(L1Catalog$Accession, 
                                function (x) length(grep(x, BestPrimerReordered$Name)) > 0)
sum(L1Catalog$WithPrimers & (!blnInReference), na.rm = T)
sum(L1Catalog$WithPrimers)
sum(L1Catalog$WithPrimers & (blnInReference), na.rm = T)
sum(blnInReference, na.rm = T)
sum(!blnInReference, na.rm = T)

# Save primer set as csv file
write.csv(as.matrix(BestPrimerReordered), 
          "D:/L1polymORF/Data/L1CatalogPrimers.csv", row.names = F)


# Count number of replicates for each name and retain only un-replicated names
# NrRep <- table(Info_isPCR$Name)
# AccNonRep <- names(NrRep)[NrRep == 1]
# Info_isPCR <- Info_isPCR[Info_isPCR$Name %in% AccNonRep, ]

# # Match the two infos
# InfoMatch <- match(Info_isPCR$Name, Info_P3$Name)
# Info_P3   <- Info_P3[InfoMatch, ]
# 
# 
# 
# 
# # Weed out invalid primers
# blnRightLength <- Info_isPCR$LengthNum <= 500 & Info_isPCR$LengthNum >= 400
# blnRightChrom  <- L1Catalog$Chromosome[AccMatch] == Info_isPCR$Chromosome
# blnLeftCovered    <- L1Catalog$start_HG38[AccMatch] > Info_isPCR$start_HG38 &
#   L1Catalog$start_HG38[AccMatch] < Info_isPCR$end_HG38 
# blnLeftCovered    <- L1Catalog$start_HG38[AccMatch] > Info_isPCR$start_HG38 &
#   L1Catalog$start_HG38[AccMatch] < Info_isPCR$end_HG38 
# blnRightCovered    <- L1Catalog$end_HG38[AccMatch] > Info_isPCR$start_HG38 &
#   L1Catalog$end_HG38[AccMatch] < Info_isPCR$end_HG38 
# blnRightPos <- blnLeftCovered | blnRightCovered
# 
# # Subset the P3Info file to get all primers that give products of the 
# # right length and right positions
# Info_P3Subset <- Info_P3[blnRightPos & blnRightLength, ]
# all(blnInReference[AccMatch][blnRightPos & blnRightLength], na.rm = T)
# 
# # Get all sequence names
# SeqNames <- sapply(as.character(Info_P3Subset$Name), function(x){
#   SplitName <- strsplit(x, "_")[[1]]
#   paste(SplitName[1:2], collapse = "_")
# })
# 
# # Select for each sequence name the best primer
# BestPrimersRef <- data.frame(stringsAsFactors = F)
# for (SeqName in unique(SeqNames)){
#   BestPrimer  <- GetBestPrimer(SeqName, Pr3Info = Info_P3Subset)
#   BestPrimersRef <- rbind(BestPrimersRef, BestPrimer)
# }
# 
# # Double check that the product length is the same as the one found by isPCR
# NameMatch <- match(BestPrimersRef$Name, Info_isPCR$Name)
# SizeDiff  <- abs(BestPrimersRef$ProductSize - Info_isPCR$LengthNum[NameMatch])
# max(SizeDiff, na.rm = T)
# sum(SizeDiff == 0, na.rm = T)
# BestPrimersRef$ProductSizeAlt <- NA
# 
# # Get all the L1s that are not in the reference
# Info_P3NotRef <- Info_P3[Info_P3$Accession %in% AccNotInRef, ]
# Info_P3NotRef$Name %in% Info_isPCR$Name
# blnDiffLength <- Info_isPCR$LengthNum <= 100 | Info_isPCR$LengthNum >= 1000
# InfoSubset <- Info_isPCR[blnDiffLength, ]
# Info_P3NotRefCandidates <- Info_P3NotRef[Info_P3NotRef$Name %in% InfoSubset$Name,]
# 
# # Get the sequence name of 
# SeqNamesNonRef <- sapply(as.character(Info_P3NotRefCandidates$Name), function(x){
#   SplitName <- strsplit(x, "_")[[1]]
#   paste(SplitName[1:2], collapse = "_")
# })
# 
# # Select for each sequence name the best primer
# BestPrimersNonRef <- data.frame(stringsAsFactors = F)
# for (SeqName in unique(SeqNamesNonRef)){
#   BestPrimer        <- GetBestPrimer(SeqName, Pr3Info = Info_P3NotRefCandidates)
#   BestPrimersNonRef <- rbind(BestPrimersNonRef, BestPrimer)
# }
# 
# # Get the alternative product size
# NameMatchNonRef <- match(BestPrimersNonRef$Name, Info_isPCR$Name) 
# BestPrimersNonRef$ProductSizeAlt <- Info_isPCR$Length[NameMatchNonRef]
# 
# # Append data set to each other
# BestPrimers <- rbind(BestPrimersRef, BestPrimersNonRef)
# write.csv(BestPrimers, "D:/L1polymORF/Data/L1CatalogPrimers.csv", row.names = F)
# 
# # Determine for each accession in L1 how often it is covered in BestPrimers
# sapply(L1Catalog$Accession, function(x) sum(BestPrimers$Accession == x, na.rm = T))

# The following script compares simulated and estimated L1 and adjusts
# the estimated L1 width to account for detection inaccuracoes

# Load packages
library(GenomicRanges)
library(MASS)
library(fields)
library(raster)

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

##########################################
#                                        #
#     Define functions                   #
#                                        #
##########################################

# Function to add columns to L1 detection dataframe
AddCols2L1Detect <- function(L1DetectDF){
  L1DetectDF$L1StartTrue    <- as.numeric(as.character(L1DetectDF$L1StartTrue))
  L1DetectDF$L1EndTrue      <- as.numeric(as.character(L1DetectDF$L1EndTrue))
  L1DetectDF$blnPass        <- L1DetectDF$EstFilter == "PASS"
  L1DetectDF$blnDetectPass       <- L1DetectDF$blnDetect & L1DetectDF$blnPass
  L1DetectDF$L1WidthFromStartEnd <- L1DetectDF$L1EndEst - L1DetectDF$L1StartEst
  L1DetectDF$L1WidthAbsDiff      <- abs(L1DetectDF$L1widthTrue - 
                                        L1DetectDF$L1widthEst)
  L1DetectDF$L1WidthStEAbsDiff   <- abs(L1DetectDF$L1widthTrue - 
                                        L1DetectDF$L1WidthFromStartEnd)
  L1DetectDF$L1StartAbsDiff      <- abs(L1DetectDF$L1StartTrue - 
                                          L1DetectDF$L1StartEst)
  L1DetectDF$L1EndAbsDiff   <- abs(L1DetectDF$L1EndTrue - L1DetectDF$L1EndEst)
  L1DetectDF$L1PosAbsDiff   <- abs(L1DetectDF$PosTrue - L1DetectDF$PosEst)
  L1DetectDF$blnFullTrue    <- L1DetectDF$L1widthTrue >= 6000
  L1DetectDF$blnFullEst     <- L1DetectDF$L1widthEst >= 6000
  L1DetectDF
}

# Function to aggregate L1Detect
AggL1Detect <- function(L1DetectDF, L1PosRange = 30){
  
  # Create an ID per L1 insertion (close insertions are considered one)
  L1DetectDF$L1ID  <- CollapseClosePos_idx(DF = L1DetectDF, 
                                           ChromCol = "Chrom", 
                                           PosCol = "PosEst", 
                                           OLRange = L1PosRange, 
                                           blnPairwise = F)
  
  # Aggregate by L1 ID
  L1DetectAgg <- AggDataFrame(L1DetectDF, 
                              GroupCol = "L1ID", 
                              MeanCols = c("L1widthEst", "L1StartEst", "L1EndEst"),
                              SumCols = c("L1GenoEst", "L1GenoTrue"),
                              MedCols = c("L1widthEst", "L1StartEst", "L1EndEst"),
                              MaxCols = "L1EndEst",
                              MinCols = "L1StartEst",
                              LengthCols = "L1widthEst",
                              Addcols = c("Chrom", "PosTrue", "L1widthTrue", 
                                          "L1StartTrue", "L1EndTrue"))
  
  # Add additional columns
  L1DetectAgg$L1Width_minmax     <- L1DetectAgg$L1EndEst_max - 
    L1DetectAgg$L1StartEst_min
  L1DetectAgg$L1Width_mix        <- L1DetectAgg$L1Width_minmax
  L1DetectAgg$L1Width_mix[L1DetectAgg$L1Width_minmax < 1000] <- 
    L1DetectAgg$L1widthEst_med[L1DetectAgg$L1Width_minmax < 1000]
  L1DetectAgg$L1WidthAbsDiff_med    <- abs(L1DetectAgg$L1widthTrue - 
                                             L1DetectAgg$L1widthEst_med)
  L1DetectAgg$L1WidthAbsDiff_minmax <- abs(L1DetectAgg$L1widthTrue - 
                                             L1DetectAgg$L1Width_minmax)
  L1DetectAgg$L1WidthAbsDiff_mix <- abs(L1DetectAgg$L1widthTrue - 
                                          L1DetectAgg$L1Width_mix)
  L1DetectAgg
  
  }  

##########################################
#                                        #
#            Load data                   #
#                                        #
##########################################

# Load data and add columns
load("D:/L1polymORF/Data/L1Simulated_MELT.RData")
L1Detect       <- AddCols2L1Detect(L1Detect)
L1Detect_Group <- AddCols2L1Detect(L1Detect_Group)

# Add numeric columns for L1 start and end
L1_1000G$L1StartNum <- as.numeric(as.character(L1_1000G$L1Start))
L1_1000G$L1EndNum   <- as.numeric(as.character(L1_1000G$L1End))

# Add columns for estimated L1 length to L1DetectAgg_Group
L1DetectAgg_Group$L1widthEst <- sapply(L1DetectAgg_Group$INFO, GetFromVcfINFO_SVLength)
L1StartEnd  <- t(sapply(VcfFile$INFO, GetFromVcfINFO_MELT_L1StartEnd))
L1DetectAgg_Group$L1StartEst <- L1StartEnd[,1]
L1DetectAgg_Group$L1EndEst   <- L1StartEnd[,2]
L1DetectAgg_Group$L1WidthAbsDiff  <- abs(L1DetectAgg_Group$L1widthTrue - 
                                           L1DetectAgg_Group$L1widthEst)

# Add frequency columns
L1DetectAgg_Group$TrueFreq <- rowSums(L1DetectAgg_Group[,as.character(SampleIDs)], na.rm = T)
idxEstGenoCol <- grep("hg19_", colnames(L1DetectAgg_Group))
for (i in idxEstGenoCol){
  L1DetectAgg_Group[,i] <- sapply(L1DetectAgg_Group[,i], GetFromVcfGeno_GenoNum)
}
L1DetectAgg_Group$EstFreq  <- rowSums(L1DetectAgg_Group[, idxEstGenoCol], na.rm = T)

# Save data with additional info
save.image("D:/L1polymORF/Data/L1Simulated_AdditionalInfo_MELT.RData")


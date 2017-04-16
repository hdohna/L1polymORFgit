# The script below reads the L1 catalog and calculates some basic numbers

##############
# Source prerequisites
##############

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

######
# Read and process L1 catalog
######

# Read in table with known L1 
L1Catalog <- read.csv("D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv", 
                        as.is = T)
L1Catalog$InRef <- L1Catalog$end_HG38 - L1Catalog$start_HG38 >= 6000

# Explore how many could not be located by author
table(is.na(L1Catalog$end_HG38), L1Catalog$Reference)

# Explore how many are in the reference genome by author
table(L1Catalog$InRef, L1Catalog$Reference)

L1Catalog[!L1Catalog$InRef & L1Catalog$Reference == "Brouha2003", ]
L1Catalog[which(L1Catalog$InRef & L1Catalog$Reference == "Beck2010"), ]

# Explore activity by author
ActivityNum <- L1Catalog$Activity
ActivityNum <- gsub("<", "", ActivityNum)
ActivityNum <- as.numeric(ActivityNum)
L1Catalog$ActivityNum <- ActivityNum

table(L1Catalog$Activity > 0, L1Catalog$Reference)

L1CatalogActive <- subset(L1Catalog, subset = ActivityNum > 0 & Allele == 1)
AlleleFreqNum <- as.numeric(L1CatalogActive$Allele_frequency)
meanFreq <- mean(AlleleFreqNum, na.rm = T)
meanFreq * nrow(L1CatalogActive)

L1CatalogSubset <- L1Catalog[L1Catalog$ActivityNum > 0,] 
aggregate(ActivityNum~Reference, data = L1Catalog, FUN = max)
aggregate(ActivityNum~Reference, data = L1CatalogSubset, FUN = min)

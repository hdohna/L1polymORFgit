# The following script reads in a version of the L1 catalog, adds columns
# and saves the result as "L1CatalogExtended.csv:

# Read in table with known L1 
L1Catalog <- read.csv("D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_16-33-31_2016.csv", as.is = T)

# Add numeric activity
L1Catalog$ActivityNum <- L1Catalog$Activity
L1Catalog$ActivityNum <- gsub("<", "", L1Catalog$ActivityNum)
L1Catalog$ActivityNum <- as.numeric(L1Catalog$ActivityNum)

# Add boolean variable indicating L1s in reference genopme
L1Catalog$blnInRef <- (L1Catalog$end_HG38 - L1Catalog$start_HG38) > 6000 

# Write out extended catalog
write.csv(L1Catalog, "D:/L1polymORF/Data/L1CatalogExtended.csv")
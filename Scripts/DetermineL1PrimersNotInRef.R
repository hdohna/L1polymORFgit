# The script below determines which primers are not in the reference

# Read in L1 catalog 
L1Catalog <- read.csv("D:/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv",
                      as.is = T)

# Indicator whether L1 is in reference
blnInReference <- (L1Catalog$end_HG38 - L1Catalog$start_HG38) > 5000
AccNotInRef  <- L1Catalog$Accession[which(!blnInReference)]
 
# Read in primer file
L1primers <- read.csv("D:/L1polymORF/Data/L1CatalogPrimers.csv",
                      as.is = T)

# Get accession numbers from primers
AccFromPrimers <- sapply(L1primers$Name, function(x) strsplit(x, "_")[[1]][1])
AccFromPrimers <- unique(AccFromPrimers)
sum(AccFromPrimers %in% AccNotInRef)
sum(!AccFromPrimers %in% AccNotInRef)

# Read file with primer results
L1_PCR_results <- read.delim("D:/L1polymORF/Data/L1_PCR_results.txt",
   as.is = T, skip = 1, sep = " ", col.names = c("Name", "Well", "Product"))

# Get accession number from file with primer results and match to catalog
L1_PCR_results$AccNr <- sapply(L1_PCR_results$Name, function(x) strsplit(x, "_")[[1]][1])
AccMatch <- match(L1_PCR_results$AccNr, L1Catalog$Accession)

# Determine which row from PCR results is in reference
L1_PCR_results$inRef <- blnInReference[AccMatch]
table(L1_PCR_results$Product, L1_PCR_results$inRef)

# Look at which coriell ID yielded a product
table(L1Catalog$Coriell_ID[AccMatch], L1_PCR_results$Product)

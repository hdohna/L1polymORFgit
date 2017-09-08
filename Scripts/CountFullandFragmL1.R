# The following script counts the full-length and fragment L1 in reference and
# 1000 genome data

# Load 1000 genome data
load('D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData')

# Read in L1 table from reference genome
L1Ref <- read.csv("D:/L1polymORF/Data/repeatsHg19_L1HS.csv")

# Number of full-length and fragment L1s in 1000 genome data
L1Counts_1000G = c(
  Full = sum(L1_1000G_reduced$InsLength > 6000, na.rm = T),
  Fragm = sum(L1_1000G_reduced$InsLength <= 6000, na.rm = T),
  LengthUnknown = sum(is.na(L1_1000G_reduced$InsLength)),
  IntactORF = NA
)

# Number of full-length and fragment L1s in reference data
L1Ref$Width <- L1Ref$genoEnd - L1Ref$genoStart
L1Counts_Ref = c(
  Full = sum(L1Ref$Width > 6000, na.rm = T),
  Fragm = sum(L1Ref$Width <= 6000, na.rm = T),
  LengthUnknown = sum(is.na(L1Ref$Width)),
  IntactORF = 90
)

# Create a data.frame with both columns and save it
L1Overview <- data.frame(Counts_1000Genomes = L1Counts_1000G,
                         Counts_ReferenceGenome = L1Counts_Ref)
write.csv(L1Overview, "D:/L1polymORF/Data/L1Overview.csv")


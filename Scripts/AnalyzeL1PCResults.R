# The following script test for each L1 in the catalog whether it is present
# in a genome using PacBio data from a capture experiment

# Source start script
#source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(csaw)
library(ShortRead)
library(Rsamtools)

# Specify flank width
Fwidth <- 50

# Specify path to PacBio bam file (reads alinged to 
# L1CatalogueWithFlank_Sat_May_07_15-15-31_2016.fa)
BamFilePath      <- "D:/L1polymORF/Data/NA12878_L1cataloguePCR.sorted.bam"
list.files("D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_ROI_catl110kb/")

# Read in table with known L1 
L1Catalog <- read.csv("D:/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv",
                      as.is = T)

# Load data with L1 location information within the catalog
load('D:/L1polymORF/Data//L1CatalogueWithFlank_Sat_May_07_15-15-31_2016_L1Locations.RData')
L1Ranges   <- GRanges(seqnames = colnames(L1StartEnd),
                             ranges = IRanges(start = L1StartEnd["Start", ], 
                                              end = L1StartEnd["End", ]))
export.bed(L1Ranges,'D:/L1polymORF/Data/L1RangesCatWithFlank_Sat_May_07_15-15-31_2016.bed')

# Determine sequence length and get all read ranges
RefSeqLengths <- sapply(L1withFlank, length)
names(L1withFlank)

# Loop through catalog elements and extract ranges
cat("Extracting reads from bam file ... \n")
ReadRanges <- lapply(1:length(RefSeqLengths), function(i){
  cat("Extracting reads for L1",names(L1withFlank)[i],"\n")
  GR <- GRanges(seqnames = names(L1withFlank)[i], 
                ranges = IRanges(start = 1, end = RefSeqLengths[i]))
  try(extractReads(BamFilePath, GR))
})
blnGRanges <- sapply(ReadRanges, function(x)is(x)[1] == "GRanges")
ReadRanges <- GRangesList(ReadRanges[blnGRanges])
ReadRanges <- unlist(ReadRanges)

# Define flank ranges for getting reads intersecting with
LeftFlankRanges   <- GRanges(seqnames = colnames(L1StartEnd),
                             ranges = IRanges(start = L1StartEnd["Start", ] - 300, 
                                              end = L1StartEnd["Start", ] - 5))
RightFlankRanges  <- GRanges(seqnames = colnames(L1StartEnd),
                             ranges = IRanges(start = L1StartEnd["End", ] + 5, 
                                              end = L1StartEnd["End", ] + 300))

# Count overlaps for left and right 
NrOverlap_Left  <- countOverlaps(LeftFlankRanges, ReadRanges, minoverlap = 1)
NrOverlap_Right <- countOverlaps(RightFlankRanges, ReadRanges, minoverlap = 1)
ReadRanges[NrOverlap_Left]

# Match L1 catalog to ranges
AccMatch <- match(as.vector(seqnames(LeftFlankRanges)), L1Catalog$Accession)
L1Catalog <- L1Catalog[AccMatch,]

# Get nor overlap for 5' and 3' end
blnNeg <- L1Catalog$Strand == "-"
NrOverlap_5P <- NrOverlap_Left
NrOverlap_5P[blnNeg] <- NrOverlap_Right[blnNeg]
NrOverlap_3P <- NrOverlap_Right
NrOverlap_3P[blnNeg] <- NrOverlap_Left[blnNeg]

hist(NrOverlap_5P, breaks = seq(0, max(NrOverlap_5P) + 10^2, 10^2))
#hist(NrOverlap_3P, breaks = 0:max(NrOverlap_3P))

# Read file with PCR results
L1_PCR_results <- read.delim("D:/L1polymORF/Data/L1_PCR_results.txt", header = F,
   as.is = T, skip = 2, sep = " ", col.names = c("Name", "Well", "Product"))

# Get accession number from file with primer results and match to catalog
L1_PCR_results$AccNr <- sapply(L1_PCR_results$Name, function(x) strsplit(x, "_")[[1]][1])
AccMatch <- match(L1_PCR_results$AccNr, colnames(L1StartEnd))
NrOverlap_5P_matched <- NrOverlap_5P[AccMatch]
NrOverlap_3P_matched <- NrOverlap_3P[AccMatch]
NrOverlap_Left_matched  <- NrOverlap_Left[AccMatch]
NrOverlap_Right_matched <- NrOverlap_Right[AccMatch]

# Get for each L1_PCR_results the number of matching reads
NcharAcc <- nchar(L1_PCR_results$AccNr)
L1_PCR_results$Side <- substr(L1_PCR_results$Name, NcharAcc + 3, NcharAcc + 3)
NrOverlap_matched <- sapply(1:nrow(L1_PCR_results), function(x){
  if (L1_PCR_results$Side[x] == "1"){
    NrOverlap_Left_matched[x]
  } else {
    NrOverlap_Right_matched[x]
  }
})
L1_PCR_results$NrOverlap <- NrOverlap_matched

# Analyze relationship between number of overlapping reads and PCR product
NrOverlap_5P_matched
table(L1_PCR_results$Product, NrOverlap_5P_matched > 0)
table(L1_PCR_results$Product, NrOverlap_3P_matched > 0)
table(L1_PCR_results$Product, NrOverlap_5P_matched > 0 | NrOverlap_3P_matched > 0)
table(L1_PCR_results$Product, NrOverlap_matched > 0)
table(L1_PCR_results$Product, NrOverlap_matched > 1000)
boxplot(NrOverlap_5P_matched ~ L1_PCR_results$Product)
boxplot(NrOverlap_3P_matched ~ L1_PCR_results$Product)
boxplot(NrOverlap_matched ~ L1_PCR_results$Product)

# Get all 

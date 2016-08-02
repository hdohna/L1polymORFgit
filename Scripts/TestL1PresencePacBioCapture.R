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
Fwidth <- 500

# Specify path to PacBio bam file
BamFilePath      <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/BZ_NA12878L1capt5-9kb_sub_reads_aln2Cat10000.sorted.bam"
BamFilePath      <- "D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_sub_reads_aln2Cat10000.sorted.bam"
BamFilePath      <- "D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_ROI_catl110kb/BZ_NA12878L1capt5-9kb_ROI_catl110kb.bam"
list.files("D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_ROI_catl110kb/")
# Read in table with known L1 
# L1Catalog <- read.csv("/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv",
#                         as.is = T)
L1Catalog <- read.csv("D:/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv",
                      as.is = T)

# Load data with L1 location information within the catalog
#load('/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/L1CatalogueWithFlank10000_Sat_May_07_15-15-31_2016_L1Locations.RData')
load('D:/L1polymORF/Data/L1CatalogueWithFlank10000_Sat_May_07_15-15-31_2016_L1Locations.RData')

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

# Get flanks and count overlaps with flanks 
# Define flank ranges for getting reads intersecting with 
LeftFlankRanges   <- resize(StartEndGR, width = 1, fix = "start")
LeftFlankRanges   <- resize(LeftFlankRanges, width = Fwidth, fix = "center")
RightFlankRanges  <- resize(StartEndGR, width = 1, fix = "end")
RightFlankRanges  <- resize(RightFlankRanges, width = Fwidth, fix = "center")


# Count overlaps for left and right 
NrOverlap_Left  <- countOverlaps(LeftFlankRanges, ReadRanges, minoverlap = Fwidth)
NrOverlap_Right <- countOverlaps(RightFlankRanges, ReadRanges, minoverlap = Fwidth)

# Match L1 catalog to ranges
AccMatch <- match(as.vector(seqnames(StartEndGR)), L1Catalog$Accession)
L1Catalog <- L1Catalog[AccMatch,]

# Get nor overlap for 5' and 3' end
blnNeg <- L1Catalog$Strand == "-"
NrOverlap_5P <- NrOverlap_Left
NrOverlap_5P[blnNeg] <- NrOverlap_Right[blnNeg]
NrOverlap_3P <- NrOverlap_Right
NrOverlap_3P[blnNeg] <- NrOverlap_Left[blnNeg]

hist(NrOverlap_5P, breaks = 0:300)
hist(NrOverlap_3P, breaks = 0:600)

# Read file with PCR results
L1_PCR_results <- read.delim("D:/L1polymORF/Data/L1_PCR_results.txt", header = F,
   as.is = T, skip = 1, sep = " ", col.names = c("Name", "Well", "Product"))

# Get accession number from file with primer results and match to catalog
L1_PCR_results$AccNr <- sapply(L1_PCR_results$Name, function(x) strsplit(x, "_")[[1]][1])
AccMatch <- match(L1_PCR_results$AccNr, as.vector(seqnames(StartEndGR)))
NrOverlap_5P_matched <- NrOverlap_5P[AccMatch]
NrOverlap_3P_matched <- NrOverlap_3P[AccMatch]

# Analyze relationship between number of overlapping reads and PCR product
NrOverlap_5P_matched
table(L1_PCR_results$Product, NrOverlap_5P_matched > 0)
table(L1_PCR_results$Product, NrOverlap_3P_matched > 0)
table(L1_PCR_results$Product, NrOverlap_5P_matched > 0 & NrOverlap_3P_matched > 0)
boxplot(NrOverlap_5P_matched ~ L1_PCR_results$Product)
boxplot(NrOverlap_3P_matched ~ L1_PCR_results$Product)

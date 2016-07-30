# The following script test for each L1 in the catalog whether it is present
# in a genome using PacBio data from a capture experiment

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(csaw)
library(ShortRead)
library(Rsamtools)

# Specify path to PacBio bam file
BamFilePath      <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/BZ_NA12878L1capt5-9kb_sub_reads_aln2Cat10000.sorted.bam"

# Read in table with known L1 
L1Catalog <- read.csv("/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv",
                        as.is = T)
# Load data with L1 location information within the catalog
load('/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/L1CatalogueWithFlank10000_Sat_May_07_15-15-31_2016_L1Locations.RData')


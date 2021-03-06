# The following script reads a bam file of reads aligned to a catalogue of
# full-length L1 and filters reads by various quality criteria

#############################################
#                                           #
#    Source packages and set parameters     #
#                                           #
#############################################

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(csaw)
library(ShortRead)
library(Rsamtools)

# Specify file paths
#BamFile         <- '/share/diskarray3/hzudohna/NA12878PacBio_aln2Catalogue2016-05-07sortedunique.sorted.bam'
BamFile         <- '/share/diskarray3/hzudohna/NA12878PacBio_aln2Catalogue2016-05-07filteredByLength.bam'
AlignListFile   <- "/home/hzudohna/L1polymORF/Data/L1CatalogueWithFlank_Sat_May_07_15-15-31_2016_L1Locations.RData"
OutFilePrefix   <- "/share/diskarray3/hzudohna/PacBio/L1PacBioLengthPosFiltered_"

# Specify parameter for border width in bp
BorderWidth  <- 50

############################
#                          #
#        Read Data         #
#                          #
############################

load(AlignListFile)

cat("Filtering bam file ...\n")
for (i in 1:ncol(L1StartEnd)){
  cat("Filtering element", i, "of", ncol(L1StartEnd), "\n")
  
  FRange <- GRanges(seqnames = colnames(L1StartEnd)[i], 
                    IRanges(start = 1, end = 16000))
  L1Border <- c((L1StartEnd[1,i] - BorderWidth):(L1StartEnd[1,i] + BorderWidth),
                (L1StartEnd[2,i] - BorderWidth):(L1StartEnd[2,i] + BorderWidth))
  PosFilter <- FilterRules(getIDs <- function(DF){!DF$pos %in% L1Border})
  paramFilter  <- ScanBamParam(which = FRange, what = scanBamWhat())
  OutFile <- paste(OutFilePrefix, colnames(L1StartEnd)[i], ".bam", sep = "")
  filterBam(BamFile, OutFile, param = paramFilter,  filter = PosFilter)
}
cat("...Done!")

# Merge filtered files
cat("Merging filtered bam files ...\n")

# Get folder, file prefix and files to be merged
OutFilePrefix_split   <- strsplit(OutFilePrefix, "/")[[1]]
OutFolder <- paste(OutFilePrefix_split[-length(OutFilePrefix_split)], 
                   collapse = "/")
FilePrefix <- OutFilePrefix_split[length(OutFilePrefix_split)]
FilesToMerge <- list.files(OutFolder, pattern = FilePrefix, full.names = T)
FilesToMerge <- FilesToMerge[-grep("bam.", FilesToMerge)]

# Merge files
MergedFile <- paste(OutFilePrefix, "Merged.bam", sep = "")
mergeBam(FilesToMerge, MergedFile, indexDestination = T)





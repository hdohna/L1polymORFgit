##############################################
#
# General description:
#
#   The following script reads a L1 catalog and creates a separate alignment for
#   all reads that overlap the junction.

# Input:
#
#     BamFile: path to file that contains mapped reads
#     L1TableFileName: path to file that contains L1HS ranges in a table

# Output:
#   
#    : ...

##############################################

######                                      
# Source packages and set parameters  
######

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(ShortRead)
library(rtracklayer)
library(Rsamtools)
library(csaw)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

# Files and folders
L1CataloguePath <- "D:/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv"
ChainFile       <- 'D:/L1polymORF/Data/hg38ToHg19.over.chain'
BamFilePath     <- 'D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_subreads_hg19.bam'

# Specify the minimum read depth to create alignment
MinReadDepth <- 5

#######                       
# Get L1 ranges                    
#######                                     

cat("Getting reference L1 ranges \n")

# Read in table with known L1 
L1Catalogue <- read.csv(L1CataloguePath, as.is = T)

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
blnInRef          <- L1Catalogue$end_HG38 - L1Catalogue$start_HG38 > 5900
L1CatalogL1Mapped <- L1Catalogue[which(blnL1Mapped & blnAllele1 & blnInRef),]

# Lift coordinates and get genomic ranges for catalog L1 on hg19
LiftOverList <- LiftoverL1Catalog(L1CatalogL1Mapped,
    ChainFilePath = ChainFile)
L1CatalogGR <- LiftOverList$GRCatalogue_hg19
L1CatalogGR <- resize(L1CatalogGR, 20000, fix = "center")

#####                                   
# Write IDs of reads in L1Ranges                        
#####

cat("Filtering bam file", BamFilePath, "\n")
R <- GRanges(seqnames = "chr1", 
             ranges = IRanges(start = 210091000, end = 210091000))
paramScan  <- ScanBamParam(which = R, what = scanBamWhat())
Reads <- scanBam(BamFilePath, param = paramScan)
Reads[[1]]$pos
Reads[[1]]$strand

# Get common start and end for all reads.
StartAll <- max(Reads[[1]]$pos)
EndAll   <- min(Reads[[1]]$pos + width(Reads[[1]]$seq))
RelStart <- StartAll - Reads[[1]]$pos + 1
RelEnd   <- RelStart + EndAll - StartAll - 1

# Get sequences in common range 
Seqs <- Reads[[1]]$seq
for (i in 1:length(Seqs)){
  if (Reads[[1]]$strand[i] == "-"){
    Seqs[i] <- reverseComplement(Seqs[i])
  }
  Seqs[i] <- as(Seqs[i][[1]][RelStart[i]:RelEnd[i]], "DNAStringSet")
}

# Get the reference sequence in common range
RefSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, names = Reads[[1]]$rname[1], 
                 start = StartAll, end = EndAll)
SeqsWithRef <- as(RefSeq, "DNAStringSet")
writeFastac(writeFasta, )


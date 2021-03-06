##############################################
#
# General description:
#
#   The following script reads in L1 catalog and looks in bam file for reads
#   aligning to L1 flanks, and saves a list of read IDs to an RData file.

# Input:
#
#  Folders and file names:
#     BamFilePath : path to bam file sith aligned reads
#     L1CatalogPath :   path to table with L1 insertions
#     L1HSConsensusPath : path to fasta file containing the L1HS consensus sequence
#     OutputFilePath:  path to file where analysis data output is saved

#  Peak calling parameters:
#     FlankWidth: width of L1 flank to get reads from 


# Output:
#   
#    RData file with a list of read IDs aligning to L1 flanks

##############################################


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
BamFilePath        <- '/share/diskarray4/MEI/NA12878/mergeBAM/NA12878.sorted.dedup.bam'
L1CatalogPath      <- '/home/hzudohna/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv'
L1HSConsensusPath  <- "/home/hzudohna/L1polymORF/Data/Homo_sapiens_L1_consensus.fa"
OutputFilePath     <- '/share/diskarray3/hzudohna/NA12878shortInsert/L1Catalog_07_15-15-31_2016_FlankIDs.RData'

# Set other parameters
FlankWidth <- 500

#############################################
#                                           #
#    Determine flank ranges                 #
#                                           #
#############################################

# Read in catalog file
L1Catalog <- read.csv(L1CatalogPath, as.is = T)

# Determine L1s not in reference
idxNonRef <- which((L1Catalog$end_HG38 - L1Catalog$start_HG38) < 5500)
L1CatalogNonRef <- L1Catalog[idxNonRef,]

# Lift coordinates to hg19
GRCatalogue_hg38  <- GRanges(seqnames = L1CatalogNonRef$Chromosome,
                             ranges = IRanges(start = pmin(L1CatalogNonRef$start_HG38,
                                                           L1CatalogNonRef$end_HG38),
                                              end = pmax(L1CatalogNonRef$start_HG38,
                                                         L1CatalogNonRef$end_HG38)),
                             strand = L1CatalogNonRef$Strand)
GRCatalogue_hg19 <- liftOver(GRCatalogue_hg38, 
                             chain = import.chain(
                               "/home/hzudohna/L1polymORF/Data/hg38ToHg19.over.chain"))
NrMapped_hg19    <- sapply(GRCatalogue_hg19, length)
idxUniqueMapped  <- which(NrMapped_hg19 == 1) 
GRCatalogue_hg19 <- unlist(GRCatalogue_hg38[idxUniqueMapped])
AccessionMapped  <- L1CatalogNonRef$Accession[idxUniqueMapped]

# Create genomic ranges per non-reference L1 flank
GRCatalogueFlank_hg19 <- flank(GRCatalogue_hg19, width = FlankWidth)
SeqNVect <- as.vector(seqnames(GRCatalogueFlank_hg19))
SeqNVect <- substr(SeqNVect, 4, nchar(SeqNVect))
GRCatalogueFlank_hg19 <- GRanges(SeqNVect,
                                 IRanges(start(GRCatalogue_hg19),
                                                  end(GRCatalogue_hg19)), 
                                 strand = strand(GRCatalogue_hg19))

#############################################
#                                           #
#    find reads and map to consensus L1     #
#                                           #
#############################################

# Loop through flank ranges, find IDs of reads mapped to flank, write out 
# unmapped reads with same ID and mapp them to consensus L1
L1FlankIDList <- lapply(1:length(GRCatalogue_hg19), function(i){
  
  # Get IDs of reads mapped to flanl
  Acc <- AccessionMapped[i]
  cat("Gettin IDs of reads mapped to L1", Acc,  "...\n")
  GR <- GRCatalogueFlank_hg19[i]
  readParams <- ScanBamParam(which = GR, what = "qname")
  IDs <- scanBam(BamFilePath, param = readParams)
  IDs <- unlist(IDs)
})
names(L1FlankIDList) <- AccessionMapped

# Save L1FlankIDList
cat("\nSaving results to", OutputFilePath)
save(list = c("L1FlankIDList", "GRCatalogueFlank_hg19"), file = OutputFilePath)

scanBamFlag(isProperPair = T)
bamFlagAsBitMatrix(scanBamFlag(isProperPair = T))

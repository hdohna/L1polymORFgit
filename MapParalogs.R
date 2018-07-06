# The following script maps paralogs

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

# Paths
ChainPath <- "D:/L1polymORF/Data/hg19.hg19.processed.chain"
OutputPath <- "D:/L1polymORF/Data/Paralogs_hg19.RData"

# Maximum liftover gap (liftover ranges that are this gap apart will be collapsed)
MaxLiftGap    <- 50

# Minimum proportion of focal gene length another gene has to overlap with to
# be counted a paralog
MinGeneOLProp <- 0.8

##########################################
#                                        #
#     Analysis                           #
#                                        #
##########################################

# Define genomic ranges of genes
GeneGR <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Get the number of matching ranges within the human genome
GeneLiftOver <- liftOver(GeneGR, 
                         chain = import.chain(con = ChainPath))

# Collapse to all
GeneLiftOver_collapsed <- lapply(GeneLiftOver, function(x) {
  CollapseGRanges(x, Width2Add = 2*MaxLiftGap)
})

# Count how many ranges a gene maps to
NrMappedRange  <- sapply(GeneLiftOver_collapsed, length)

# Count the number of ranges a gene maps to that overlap with another
# gene by a length that is at least the proportion MinGeneOLProp of the 
# focal gene length
NrMappedGene   <- sapply(1:length(GeneLiftOver_collapsed), function(i) {
  x <- GeneLiftOver_collapsed[[i]]
  OLWidth <- floor(MinGeneOLProp * width(GeneGR[i]))
  sum(overlapsAny(x, GeneGR, ignore.strand = T, minoverlap = OLWidth))
})
hist(NrMappedGene)

# Save objects
save(list = c("GeneLiftOver_collapsed", "NrMappedGene", "NrMappedRange",
              "MaxLiftGap", "MinGeneOLProp"),
     file = OutputPath)

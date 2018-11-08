# The script below identifies deletions that are likely L1 polymorphis,s

# Load packages
library(GenomicRanges)

##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

# Specify file paths
DataPath        <- 'D:/L1polymORF/Data/'
L1GRPath        <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'
L1RefRangePath  <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
ChrLPath        <- 'D:/L1polymORF/Data/ChromLengthsHg19.Rdata'
VariantPath     <- 'D:/L1polymORF/Data/ALL_chr1_ME_Deletions.vcf'

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("\n\nLoading and processing data ...")

# Load previously generated objects
load(L1GRPath)
load(ChrLPath)
load(L1RefRangePath)

# Load table of variants
VcfFile <- read.table(VariantPath, col.names = colnames(L1_1000G)[1:2513], as.is = T)
VcfFile$ChromNameLong <- paste("chr", VcfFile$CHROM, sep = "")

# Extract the end of the deletion from genomic ranges
VcfFile$DelEnd <- sapply(VcfFile$INFO, function(x) {
  InfoSplit = strsplit(x, ";")[[1]]
  EndPart <- InfoSplit[substr(InfoSplit, 1, 4) == "END="]
  as.numeric(substr(EndPart, 5, nchar(EndPart)))
})

# Create genomic range for vcf file
VcfGR <- makeGRangesFromDataFrame(df = VcfFile, 
                                  seqnames.field = "ChromNameLong", start.field = "POS",
                                  end.field = "DelEnd")
VcfGR <- VcfGR[width(VcfGR) <= 7000]

# Find overlap between L1 ranges and deletions
L1_1000G_GR_hg19_Plus200 <- GRanges(seqnames = seqnames(L1_1000G_GR_hg19),
                                    IRanges(start = start(L1_1000G_GR_hg19) - 200,
                                            end = end(L1_1000G_GR_hg19) + 200))
L1DelOL <- findOverlaps(L1_1000G_GR_hg19_Plus200, VcfGR, ignore.strand = T)
L1DelOL

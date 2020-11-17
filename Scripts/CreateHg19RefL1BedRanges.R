# The following script creates a bed file for the non-polymorphic L1 ranges for the reference
# genome hg19 

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(rtracklayer)

# Load L1 reference ranges
load('D:/OneDrive - American University of Beirut/L1polymORF/Data/L1RefRanges_hg19.Rdata')
width(L1GRanges)

# Read in vcf file with MELT deletion calls
MEDelCall <- ReadVCF("D:/OneDrive - American University of Beirut/L1polymORF/Data/DEL.final_comp.vcf")
MEDelCall$chromosome <- paste("chr", MEDelCall$X.CHROM, sep = "")
MEDelCall$end <- sapply(MEDelCall$INFO, function(x){
  Split1 <- strsplit(x, ";")[[1]]
  EndPart <- grep("END=", Split1, value = T)
  EndChar <- strsplit(EndPart, "END=")[[1]][2]
  as.numeric(EndChar)
})
MEDel_GR  <- makeGRangesFromDataFrame(df = MEDelCall,
                                      start.field = "POS",
                                      end.field = "end")

# Get reference L1HS that are not polymorphic
blnL1Del        <- overlapsAny(L1GRanges, MEDel_GR)
L1GRangesNoPoly <- L1GRanges[!blnL1Del]
ChrNames <- as.vector(seqnames(L1GRangesNoPoly))
L1GRangesNoPoly <- GRanges(seqnames = substr(ChrNames, 4, nchar(ChrNames)),
                           ranges = IRanges(start = start(L1GRangesNoPoly),
                                            end = end(L1GRangesNoPoly)))
L1GRangesNoPolyFull <-  L1GRangesNoPoly[width(L1GRangesNoPoly) >= 6000]
export.bed(L1GRangesNoPoly, "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1GRangesNoPoly.bed")
export.bed(L1GRangesNoPolyFull, "D:/OneDrive - American University of Beirut/L1polymORF/Data/L1GRangesNoPolyFull.bed")

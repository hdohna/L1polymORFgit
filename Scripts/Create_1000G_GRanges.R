# The script below subsets creates GRanges for all 
# variants in the 1000 Genome data
# Data were downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')
library(GenomicRanges)

# Specify parameters
NrInfoCols <- 9
DataFolder <- "/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/"
Chrom      <- "chr1"
NrRowsPerRead <- 5*10^4

# Specify vcf file and output file based on chromosome
VcfFilePattern  <- paste("ALL.", Chrom, ".phase3", sep = "")
VcfFile <- list.files(DataFolder, pattern = VcfFilePattern, full.names = T)
VcfFile <- VcfFile[-grep("gz.tbi", VcfFile)]
OutFile <- paste(DataFolder, Chrom, "_AllVariantGRanges.RData", sep = "")

# Loop over file names, read file and append to existing
cat("Reading vcf files by chromosome\n")
VarTable <- read.table(VcfFile, header = F, skip = 1, nrows = 10)

# Initialize objects
RowsReadIn   <- 0
blnRows2Read <- T
VariantGR     <- GRanges()

while(blnRows2Read){
  
  cat("Reading in rows from row", RowsReadIn + 1, "to ...")
  # Read in parts of the table with transcript expression data
  VarTable <- read.table(VcfFile, skip = RowsReadIn + (RowsReadIn == 0), 
                         skipNul = T,
                         nrows = NrRowsPerRead, as.is = T, header = F)
  cat(RowsReadIn + nrow(VarTable), "\n")
  cat("done! \n")
  
  # Append rows to growing genomic ranges object
  NewVariantGR <- GRanges(seqnames = Chrom, IRanges(start = VarTable$V2,
                                                    end = VarTable$V2),
                          mcols = VarTable$V3)
  VariantGR <- append(VariantGR, NewVariantGR)
  cat("Kept", length(VariantGR), "genomic ranges \n")
}

# Save genomic ranges
cat("Saving genomic ranges to", OutFile, "\n")
save(list = c("VariantGR", "Chrom"), file = OutFile)

cat("done!\n")




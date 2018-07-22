# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(GenomicRanges)

# Process Gencode genes
cat("Processing Gencode genes ...")
GencodeGenes   <- read.delim("D:/L1polymORF/Data/GENCODE_V19_genes.txt")
GeneGR_Gencode <- makeGRangesFromDataFrame(
  GencodeGenes[!duplicated(GencodeGenes$name2),
               c("chrom", "strand", "txStart", "txEnd", "name2")],
  keep.extra.columns = T, 
  start.field = "txStart",
  end.field   = "txEnd")
#GeneGR_Gencode <- UniqueGRanges(GeneGR_Gencode, Group = GencodeGenes$name2)
cat("done!\n")

# Process NCBI refseq genes
cat("Processing NCBIRefSeq genes ...")
NCBIRefSeqGenes <- read.delim("D:/L1polymORF/Data/NCBI_RefSeq_gene.txt")
GeneGR_NCBIRefSeq <- makeGRangesFromDataFrame(
  NCBIRefSeqGenes[!duplicated(NCBIRefSeqGenes$name2),
               c("chrom", "strand", "txStart", "txEnd", "name2")],
  keep.extra.columns = T, 
  start.field = "txStart",
  end.field   = "txEnd")
cat("done!\n")

# Process NCBI refseq genes
cat("Processing UCSC genes ...")
UCSCGenes <- read.delim("D:/L1polymORF/Data/UCSC_genes.txt")
UCSCGenes <- UCSCGenes[UCSCGenes$proteinID != "", ]
GeneGR_UCSC <- makeGRangesFromDataFrame(
  UCSCGenes[!duplicated(UCSCGenes$proteinID),
                  c("chrom", "strand", "txStart", "txEnd", "proteinID")],
  keep.extra.columns = T, 
  start.field = "txStart",
  end.field   = "txEnd")
cat("done!\n")

# Save GRanges
cat("Saving genomic ranges ...")
save(list = c("GeneGR_Gencode", "GeneGR_NCBIRefSeq",
              "GeneGR_UCSC"), file = "D:/L1polymORF/Data/Gene_GR.RData")
cat("done!\n")

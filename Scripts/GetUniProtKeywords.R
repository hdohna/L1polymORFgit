# The following script gets uniprot keywords for all genes and saves them as a csv file

# Load packages
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(UniProt.ws)

# Create genomic ranges of all genes
GeneGR <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Create a look-up table for uniprot ids
GeneLookUp <- select(org.Hs.eg.db, keys = GeneGR@elementMetadata@listData$gene_id,
       columns = "UNIPROT", keytype ="ENTREZID")
UPunique <- unique(GeneLookUp$UNIPROT)

# Create a uniprot database 
UP <- UniProt.ws(taxId = 9606)
keytypes(UP)

# Open connection to file that collects UniProt keywords and initialize
FileCon <- file("D:/L1polymORF/Data/UniProtKeywords.csv")
open(FileCon)
NrRecordsPerIter <- 500
idxStart <- seq(1, length(UPunique), NrRecordsPerIter)

# Loop and get keywords and save write them out to file
for (Start in idxStart){
  cat("Processing chunk", which(Start == idxStart), "out of", length(idxStart), "\n")
  End <- min(length(UPunique), Start + NrRecordsPerIter - 1)
  UniProtKeywords <- select(UP, keys = UPunique[Start:End],
                            columns = "KEYWORDS", keytype = "UNIPROTKB")
  UniPKeywords <- paste(UniProtKeywords$UNIPROTKB, UniProtKeywords$KEYWORDS, sep = ", ")
  writeLines(UniPKeywords, con = FileCon)
}

close(FileCon)
cat("Done!")
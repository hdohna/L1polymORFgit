##############################################
#
# General description:
#
#   The following script reads tables with accession numbers of clonal
#   sequences containing L1, downloads the sequences, writes them out
#   as fasta file and creates an index

# Input:
#
#    Brouha2003Table: table with accession numbers collected by Brouha et al 2003
#    Beck2010Table: table with accession numbers collected by Beck et al 2010
#   

# Output:
#   
#    L1CloneSeq.fas: fasta file with all clonal sequences containing L1

##############################################


# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(seqinr)

# Read in tables with clone accession numbers
Brouha2003Table <- read.csv("D:/L1polymORF/Data/L1Brouha2003.csv", as.is = T)
Beck2010Table   <- read.csv("D:/L1polymORF/Data/Beck2010_mergedTable_withEmptySite.csv",
                          as.is = T)

# Create a vector of accession numbers
AccNrs <- c(Brouha2003Table$Accession, Beck2010Table$Accession)

# Download sequences and put them in a list 
cat("***** Downloading sequences from Genbank ... \n")
choosebank("genbank")
AccNr <- AccNrs[1]
Seqs <- lapply(AccNrs, function(AccNr){
  cat("Processing", AccNr, "\n")
  x <- query(listname = "L1", paste("AC=", AccNr, sep = ""))
  Seq <- getSequence(x$req)
  toupper(Seq[[1]])
})
closebank()
cat("Done!  *****\n\n")

# Change letters to upper case and write out as fasta file
cat("***** Writing sequences to D:/L1polymORF/Data/L1CloneSeq.fas... \n")
write.fasta(Seqs, AccNrs, "D:/L1polymORF/Data/L1CloneSeq.fas")



# Load packages
library(seqinr)
library(ape)

# Read in table of hot L1 (obtained from Brouha et al.2003 PNAS)
L1Hot_raw <- read.delim('D:/L1polymORF/Data/L1HotTable_raw.txt', sep = " ")

# Find duplicated column names 
CNames        <- colnames(L1Hot_raw)
ColNameSuffix <- substr(CNames, nchar(CNames) - 1, nchar(CNames))
blnDuplCols   <- ColNameSuffix == ".1"

# Separate table according to duplicated column names and append duplucated part
L1HotDupl <- L1Hot_raw[,blnDuplCols]
colnames(L1HotDupl) <- CNames[!blnDuplCols]
L1HotTable    <- rbind(L1Hot_raw[,!blnDuplCols],L1HotDupl)

# Download sequences 
choosebank("genbank")

L1SeqHot <- lapply(L1HotTable$Accession_no., function(AccNr){
  x <- query(listname = "L1", paste("AC=", AccNr, sep = ""))
  getSequence(x$req)[[1]]
})
sapply(L1SeqHot, length)
closebank()

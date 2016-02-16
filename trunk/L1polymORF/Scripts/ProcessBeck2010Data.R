##############################################
#
# General description:
#
#   The following script reads a tables from Beck et al (2010, Cell) and 
#   assembles a dataset that contains L1 activity, genomic position and 
#   Coriell ID of carrier

# Input:
#
#    L1accession_Beck2010_table5.txt: table with fosmid and accession number
#    L1activity_Beck2010_table1.txt: table with fosmid number, chromosome, 
#       activity and L1 ID.
#    CloneID_Acc_Kidd_etal_2008.csv: table with clone name (can be linked to 
#       clone ID in CloneCoriellTable_Beck2010.txt) and GenBank accession 
#       number
#    CloneCoriellTable_Beck2010.txt: table with clone ID and CoriellID
#   

# Output:
#   
#    D:/L1polymORF/Data/Beck2010_mergedTable.csv: csv file with info merged

##############################################

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(seqinr)
library(ape)

# Read in tables
L1acc_raw    <- read.delim('D:/L1polymORF/Data/L1accession_Beck2010_table5.txt', 
                        sep = " ", as.is = T)
L1act_raw    <- read.delim('D:/L1polymORF/Data/L1activity_Beck2010_table1.txt', 
                        sep = " ", as.is = T)
CloneID_Acc  <- read.csv('D:/L1polymORF/Data/CloneID_Acc_Kidd_etal_2008.csv', 
                        as.is = T)
CloneCoriell <- read.delim('D:/L1polymORF/Data/CloneCoriellTable_Beck2010.txt', 
                          sep = " ", as.is = T)
list.files('D:/L1polymORF/Data/')

# Append duplicated columns 
L1accTable <- AppendDuplColumns(L1acc_raw)
L1actTable <- AppendDuplColumns(L1act_raw)
L1actTable <- L1actTable[-70, ]

# Create a column with short accession number
L1accTable$Accession <- sapply(L1accTable$AccessionNumber, 
                               function(x)  strsplit(x, "\\.")[[1]][1])

# Rename column "FosmidNumber" by "L1_ID"
colnames(L1accTable)[colnames(L1accTable) == "FosmidNumber"] <- "L1_ID"
MergedTable1 <- merge(L1accTable, L1actTable)

# Get the library number and merge with Colnecoriell table
MergedTable1$Library_Nr <- substr(MergedTable1$L1_ID, 1, 1)
MergedTable2 <- merge(MergedTable1, CloneCoriell, all.x = T)

# Save merged Table
write.csv(MergedTable2, "D:/L1polymORF/Data/Beck2010_mergedTable.csv")

# # Create a clone Library_ID column in CloneID_Acc
# CloneNameLong <- CloneID_Acc$Alternative.Name
# CloneNameLong <- gsub("-", "_", CloneNameLong)
# CloneName     <- sapply(CloneNameLong, function(x) strsplit(x, "_")[[1]][1])
# CloneID_Acc$Library_ID <- CloneName
# MergedTable2 <- merge(CloneID_Acc, CloneCoriell, all = T)
# 
# MergedTable3 <- merge(MergedTable1, MergedTable2)

##############################################
#
# General description:
#
#   The following script reads in retortansposn data of the euL1db data base
#   see http://eul1db.unice.fr/

# Input:
#
#    Family_eul1db: data on families
#    Individuals_eul1db: data on individuals
#    Methods_eul1db: data on methods
#    Study_eul1db: data on studies
#    Samples_eul1db: data on samples
#    MRIP_eul1db: data on methods
#    SRIP_eul1db: data on methods

# Output:
#   
#    : ...

##############################################

# Source start script
#source('D:/HumanGenome/Scripts/_Start_HumanGenome.r')

# Read in data
Families    <- read.delim("D:/HumanGenome/Data/Family_eul1db", skip = 5)
Individuals <- read.delim("D:/HumanGenome/Data/Individuals_eul1db", skip = 5)
Methods     <- read.delim("D:/HumanGenome/Data/Methods_eul1db", skip = 5)
Studies     <- read.delim("D:/HumanGenome/Data/Study_eul1db", skip = 5)
Samples     <- read.delim("D:/HumanGenome/Data/Samples_eul1db", skip = 5)
MRIP        <- read.delim("D:/HumanGenome/Data/MRIP_eul1db", skip = 5)
#SRIP        <- read.delim("D:/HumanGenome/Data/SRIP_eul1db", skip = 5)
#Reference        <- read.delim("Data/ReferenceL1HS_eul1db", skip = 5)

grep("Beck", MRIP$studies, v)

# Check study with many samples p
Samples[Samples$Sample_name == "NA12878",]
ID1878 <- Samples$Individual_id[Samples$Sample_name == "NA12878"][1]
grep(ID1878, MRIP$samples)
# Get all rows samples from individual 447
Samples447 <- which(Samples$Individual_id == 447)

# Get all SRIPs found in indvivdual 447
SRIP447 <- SRIP[SRIP$sampleID %in% Samples447,]
SRIP447 <- SRIP447[!duplicated(SRIP447$mrip),]

# determine for each sample which 
AllMrips <- SRIP447$mrip
MutMat <- sapply(Samples447, function(x) x == SRIP447$sampleID)
rownames(MutMat) <-  SRIP447$mrip
colnames(MutMat) <- Samples447
lapply(2:14, function(j) table(MutMat[,1], MutMat[,j]))


unique(SRIP$chromosome)
SRIP[(SRIP$chromosome == 1 & SRIP$g_start > 121485100 & SRIP$g_start < 121485300),]
SRIP[(SRIP$chromosome == 2 & SRIP$g_start > 92325000 & SRIP$g_start < 92325500),]

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


# Check study with many samples p
Samples[Samples$Sample_name == "NA12878",]
ID1878 <- Samples$Individual_id[Samples$Sample_name == "NA12878"][1]

idxMRIP_NA12878 <- grep(ID1878, MRIP$samples)
MRIP[idxMRIP_NA12878, c("chromosome", "start", "stop")]

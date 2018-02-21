# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Read in table with 1000 genome files
VcfFileTable <- read.table("/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/vcfFiles1000G", as.is = T)
VcfFiles <- VcfFileTable[ ,ncol(VcfFileTable)]

# Create an rsync command
i <- 1
SourcePath <- paste("rsync://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
                    VcfFiles[i], sep = "")
DestPath <- paste("./srv/gsfs0/projects/levinson/hzudohna/1000Genomes/",
                  VcfFiles[i], sep = "")
RSyncCMD <- paste("rsync -a -P", SourcePath, DestPath)

# Create and call a script file
CreateAndCallqsubScript(file = "/home/rsyncExample", 
                        qsubCommandLines = RSyncCMD)
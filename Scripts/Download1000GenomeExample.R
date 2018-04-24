# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Read in table with 1000 genome files
VcfFileTable <- read.table("D:/L1polymORF/Data/vcfFiles1000G", as.is = T)
VcfFileTable <- read.table("/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/vcfFiles1000G", as.is = T)
VcfFiles <- VcfFileTable[ ,ncol(VcfFileTable)]

# Create an rsync command
i <- 2
SourcePath <- paste("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
                    VcfFiles[i], sep = "")
DestPath <- paste("./srv/gsfs0/projects/levinson/hzudohna/1000Genomes/",
                  VcfFiles[i], sep = "")
#RSyncCMD <- paste("rsync -a -P", SourcePath, DestPath)
lftpCMD <- c("lftp",
             paste("open", SourcePath),
             paste("get", VcfFiles[i], "-o", DestPath),
             "close",
             "exit")

# Create and call a script file
CreateAndCallqsubScript(file = "/home/hzudohna/lftpExample", 
                        qsubCommandLines = lftpCMD, blnWait = T)
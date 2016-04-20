# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Source library
library(RCurl)

# URL of 1000 genome ftp site
DownloadAllFiles(DataFTP = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/',
                 FileDestination = '/share/diskarray3/hzudohna/1000Genomes')

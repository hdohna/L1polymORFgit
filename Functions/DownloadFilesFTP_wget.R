# General description:
#    The following function downloads files from an http path

# Arguments:
#   FileListPath: Path to where the list with file names are saved to.
#   FileURL : url where files are
#   BarWidth: Width of bar
#   Color: bar colors
# Output:
#   Bars added to plot

DownloadFilesFTP_wget <- function(FileURL, OutputDir = NULL, FileListPath = "FileList", 
                              FileExt = ".vcf", ScriptName = "DownloadScript",
                              RunTime = '12:00:00',
                              Mem = '50G') {
  ScriptName = "DownloadScript"
  RunTime = '12:00:00'
  Mem = '50G'
  
  source('/home/hb54/L1polymORFgit/Scripts/_Start_L1polymORF_AUB.R')
  OutputDir <- "/home/hb54/1000GenomeData"
  OutputDir <- "/scratch/hb54/dbSNP"
  OutputDir <- "/scratch/hb54/1000G_high_coverage"
  
  FileURL <- "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/"
  FileURL <- "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/"
  FileURL <- "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/annotated/"
  FileListPath <- "1000GenomeFiles"
  FileListPath <- "1000GenomeFiles_high_cover"
#  FileListPath <- "dbSNPFiles"
  
  CurlCmd <- paste("curl", FileURL, ">", FileListPath)
  system(CurlCmd)

  # Read in table with file names and write out a file with urls 
  # FileTable <- read.table(FileListPath)
  # writeLines(paste0(FileURL, FileTable[,ncol(FileTable)]), FileListPath)
  FileLines <- readLines(FileListPath)
  FileNames <- sapply(FileLines, function(x) strsplit(x, "\"")[[1]][2])
  FileNames <- FileNames[!is.na(FileNames)]
  FileNames
  grep(".vcf.gz", FileLines, value = T)
  writeLines(grep(".vcf.gz", FileLines, value = T), FileListPath)
  
  # Create wget command and run a slurm script with get command
  WgetCmd <- c(paste("cd", OutputDir), paste("wget -i", paste0("/home/hb54/", FileListPath)))
  CreateAndCallSlurmScript(file = ScriptName, 
                           scriptName = ScriptName,
                           SlurmCommandLines = WgetCmd,
                           RunTime = RunTime,
                           Mem = Mem)
  
} # End of function

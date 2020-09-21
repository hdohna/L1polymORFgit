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
  
  OutputDir <- "/home/hb54/1000GenomeData"
  FileURL <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
  FileListPath <- "1000GenomeFiles"
  
  CurlCmd <- paste("curl", FileURL, ">", FileListPath)
  system(CurlCmd)

  # Read in table with file names and write out a file with urls 
  FileTable <- read.table(FileListPath)
  writeLines(paste0(FileURL, FileTable[,ncol(FileTable)]), FileListPath)

  # Create wget command and run a slurm script with get command
  WgetCmd <- c(paste("cd", OutputDir), paste("wget -i", FileListPath))
  CreateAndCallSlurmScript(file = ScriptName, 
                           scriptName = ScriptName,
                           SlurmCommandLines = WgetCmd,
                           RunTime = RunTime,
                           Mem = Mem)
  
} # End of function

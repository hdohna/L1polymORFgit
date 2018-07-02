# The script below reads in output produced by the selscan software 
# (https://github.com/szpiech/selscan/) for the EHH option

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Determine data folder
DataFolder <- "/labs/dflev/hzudohna/1000Genomes/"
ResultPath <- "/labs/dflev/hzudohna/1000Genomes/iHHScores.RData"

# Get file names of EHH output
EHHOutFiles <- list.files(DataFolder, pattern = ".ehh.", full.names = T)
for (rmPattern in c(".colormap", ".log")){
  EHHOutFiles <- EHHOutFiles[-grep(rmPattern, EHHOutFiles)]
}

# Loop over files with EHH output and calculate iHH scores
iHHScores <- sapply(1:length(EHHOutFiles), function(i){
  
  cat("Processing result", i, "of", length(EHHOutFiles), "\n")
  # Read in the table
  FilePath <- EHHOutFiles[i]
  EHHTable <- try(read.delim(FilePath, header = F))
  L1ID     <- strsplit(FilePath, "\\.")[[1]][3]
  
  if (is(EHHTable) == "try-error"){
    iHH <- NA
  } else {
    # Get difference in genomic position and EHH values
    DiffGen  <- EHHTable$V1[-1] - EHHTable$V1[-nrow(EHHTable)] 
    DiffEHH1 <- EHHTable$V4[-1] + EHHTable$V4[-nrow(EHHTable)] 
    DiffEHH0 <- EHHTable$V5[-1] + EHHTable$V5[-nrow(EHHTable)] 
    
    # Calculate log ratio of iHH
    iHH1 <- 0.5 * DiffGen %*% DiffEHH1  
    iHH0 <- 0.5 * DiffGen %*% DiffEHH0  
    iHH  <- log(iHH1 / iHH0)
  }
  c(L1ID = L1ID, iHH = iHH)
})

# Turn matrix into data.frame
iHHScores     <- as.data.frame(t(iHHScores), stringsAsFactors = F)
iHHScores$iHH <- as.numeric(iHHScores$iHH)

# Save image 
save(ResultPath)
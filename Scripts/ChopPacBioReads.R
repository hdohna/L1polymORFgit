##############################################
#
# General description:
#
#     The following script reads reads from a PacBio fastq file and chops
#     each read in small fragments (~500bp) so that they can be aligned with
#     standard methods

# Input:
#
#     InFastQfilePath: character string indicating path to the 
#           fastq file to be read in
#     OutFilePath: character string indicating path to the 
#           fastq file to be read in
#     NrReadsPerIter: integer value indicating how many reads from the fastq
#          file should be read in per iteration
#     NewReadLength: integer value indicating nr bp per newly created read length


# Output:
#   
#    : ...

# Comments:
#    This function requires packages Rsamtools and ShortRead

##############################################

library(ShortRead)

  cat("****************************************************\n")
  cat("**                                                **\n")
  cat("**    Running script ChopPacBioReads ...     **\n")
  cat("**                                                **\n")
  cat("****************************************************\n\n")

# Set path to fastq file to be read in
InFastQfilePath <- "/share/diskarray3/hzudohna/NA12878Pacbio.fastq"  
OutFastQfilePath <- "/share/diskarray3/hzudohna/NA12878Pacbio_chopped.fastq"  
LinesPerScan <- 100000
NewReadLength <- 500

# Open one connection per fastq file and intitialize counter variables
OF <- file(OutFastQfilePath, "w")
NrReadsRead  <- 1
ReadCounter  <- 0
WriteCounter <- 0
i <- 1
# Loop through fastq file and append to little range-specific fastq files
while (NrReadsRead > 0){
  ScannedLines <- scan(InFastQfilePath, skip = (i - 1) * LinesPerScan , 
                       nlines = LinesPerScan, 
                       what = 'character') 
  
  NrReadsRead  <- length(ScannedLines) / 2
  ReadCounter  <- ReadCounter + NrReadsRead
  cat("Total of", ReadCounter, " reads read in \n")
  idxNames     <- grep(">", ScannedLines)
  if (any(idxNames %% 2 != 1)){
    stop("Reads do not start with names!")
  }
  # system.time(
  #   ReadList <- lapply(idxNames, function(i) {
  #   Read <- ScannedLines[i + 1]
  #   StartVals <- seq(1, nchar(Read), NewReadLength)
  #   ChoppedReads <- sapply(StartVals, function(x) {
  #     substr(Read, StartVals, StartVals + NewReadLength - 1)
  #   })
  #   as.vector(rbind(ScannedLines[i], ChoppedReads))
  #   })
  # )
  StartVals <- seq(1, max(nchar(ScannedLines)), NewReadLength)
  ChoppedReads <- lapply(StartVals, function(x){
    Reads <- substr(ScannedLines[idxNames + 1], x, x + NewReadLength - 1)
    ReadsWithNames <- rbind(ScannedLines[idxNames], Reads)
    ReadsWithNames[,nchar(Reads) > 0]
  })
  ChoppedReads <- unlist(ChoppedReads)
  writeLines(ChoppedReads, OF)
  WriteCounter <- WriteCounter + length(ChoppedReads) / 2
  cat("Total of", WriteCounter, " reads written out \n")
}
close(OF)

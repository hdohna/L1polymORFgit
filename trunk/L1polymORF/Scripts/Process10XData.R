# The following script processes 10X data into two fastq files

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(csaw)
library(ShortRead)
library(Rsamtools)

# Set parameters
OutFastQfilePath1 <- "/share/diskarray3/hzudohna/10XData/NA1281878_capture10X_R1.fastq"  
OutFastQfilePath2 <- "/share/diskarray3/hzudohna/10XData/NA1281878_capture10X_R2.fastq"  
LinesPerScan  <- 10^6
if(LinesPerScan %% 4 > 0){
  stop("Number of lines per scan has to be divisible by 4!")
}

# specify path to files with input data
Files <-   c("read-RA_si-TTTCATGA_lane-004-chunk-002.fastq",
             "read-RA_si-TTTCATGA_lane-003-chunk-003.fastq",
             "read-RA_si-TTTCATGA_lane-002-chunk-000.fastq",
             "read-RA_si-TTTCATGA_lane-001-chunk-001.fastq",
             "read-RA_si-GAAGGAAC_lane-004-chunk-002.fastq",
             "read-RA_si-GAAGGAAC_lane-003-chunk-003.fastq",
             "read-RA_si-GAAGGAAC_lane-002-chunk-000.fastq",
             "read-RA_si-GAAGGAAC_lane-001-chunk-001.fastq",
             "read-RA_si-CGCATGTG_lane-004-chunk-002.fastq",
             "read-RA_si-CGCATGTG_lane-003-chunk-003.fastq",
             "read-RA_si-CGCATGTG_lane-002-chunk-000.fastq",
             "read-RA_si-CGCATGTG_lane-001-chunk-001.fastq",
             "read-RA_si-ACGTCCCT_lane-004-chunk-002.fastq",
             "read-RA_si-ACGTCCCT_lane-003-chunk-003.fastq",
             "read-RA_si-ACGTCCCT_lane-002-chunk-000.fastq",
             "read-RA_si-ACGTCCCT_lane-001-chunk-001.fastq")


# Create paths to all fastq files
FilePaths <- paste("/share/diskarray2/L1HS/10X_LINE1capture/H7VK5AFXX/BCL_PROCESSOR_CS/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/",
                   Files, sep = "")

# Open one connection per fastq file and intitialize counter variables
OF1 <- file(OutFastQfilePath1, "w")
OF2 <- file(OutFastQfilePath2, "w")
NrReadsRead  <- 1
ReadCounter  <- 0
WriteCounter <- 0
i <- 1

# Loop through files and save reads to fasta files for read 1 and 2
InFastQfilePath <- "/share/diskarray2/L1HS/10X_LINE1capture/H7VK5AFXX/outs/fastq_path/read-RA_si-TTTCATGA_lane-001-chunk-001.fastq"  
for (InFastQfilePath in FilePaths){
  
  # status message for current file
  cat("******    Processing file", InFastQfilePath, "*********\n")
  
  #  Loop through fastq file and append to little range-specific fastq files
  while (NrReadsRead > 0){
     ScannedLines <- scan(InFastQfilePath, skip = (i - 1) * LinesPerScan , 
                       nlines = LinesPerScan, 
                       what = 'character', sep = '\n')
     i <- i + 1
     idxRead1 <- grep("1:N:0:0", ScannedLines)
     idxRead2 <- grep("3:N:0:0", ScannedLines)
     NrReadsRead  <- length(ScannedLines) / 4
     ReadCounter  <- ReadCounter + NrReadsRead
     cat("Total of", ReadCounter, " reads read in \n")
     idxNames     <- which(substr(ScannedLines, 1, 1) == "@")
     idxSeqs      <- idxNames + 1
     idxPlus      <- which(substr(ScannedLines, 1, 1) == "+")
     if (any(idxPlus != (idxNames + 2))) {
       stop("Lines not in proper order!")
     }
     if (any(idxNames %% 2 != 1)){
       stop("Reads do not start with names!")
     }
  
     # Add barcodes to plus lines
     NewPlus <- paste("+BX:", substr(ScannedLines[idxSeqs], 1, 16), sep = "")
     ScannedLines[idxPlus] <- NewPlus
  
     # Remove barcodes from sequences
     ScannedLines[idxSeqs] <- substr(ScannedLines[idxSeqs], 17, 
                                  nchar(ScannedLines[idxSeqs]))
  
     # Get indices for read 1 and read 2 and lines belonging to them
     idx1 <- as.vector(rbind(idxRead1, idxRead1 + 1, idxRead1 + 2, idxRead1 + 3))
     idx2 <- as.vector(rbind(idxRead2, idxRead2 + 1, idxRead2 + 2, idxRead2 + 3))
     writeLines(ScannedLines[idx1], OF1)
     writeLines(ScannedLines[idx2], OF2)
     WriteCounter <- WriteCounter + length(idxRead1)  + length(idxRead1)
     cat("Total of", WriteCounter, " reads written out \n")
   }

}

# Close file connections
close(OF1)
close(OF2)

                   
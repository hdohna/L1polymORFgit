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

# Loop through files and save reads to fasta files for read 1 and 2
InFastQfilePath <- "/share/diskarray2/L1HS/10X_LINE1capture/H7VK5AFXX/outs/fastq_path/read-RA_si-TTTCATGA_lane-001-chunk-001.fastq"  
for (InFastQfilePath in FilePaths){
  
  # status message for current file
  cat("******    Processing file", InFastQfilePath, "*********\n")
  i <- 1
  ScannedLines <- rep(1, LinesPerScan)
  
  #  Loop through fastq file and append to little range-specific fastq files
  while (length(ScannedLines)  == LinesPerScan){
     ScannedLines <- scan(InFastQfilePath, skip = (i - 1) * LinesPerScan , 
                       nlines = LinesPerScan, 
                       what = 'character', sep = '\n')
     i <- i + 1
     idxRead1 <- grep("1:N:0:0", ScannedLines)
     idxRead2 <- grep("3:N:0:0", ScannedLines)
     
     # Check that read 1 and 2 have the same names
     RN1 <- substr(ScannedLines[idxRead1], 33, 42)
     RN2 <- substr(ScannedLines[idxRead2], 33, 42)
     if (! all(RN1 == RN2)){
       stop("Names for read 1 and 2 are not consistent!")
     }
     
     # Update read counter
     NrReadsRead  <- length(ScannedLines) / 4
     ReadCounter  <- ReadCounter + NrReadsRead
     cat("Total of", ReadCounter, " reads read in \n")
 
     # Add barcodes to plus lines
     idxPlus  <- which(substr(ScannedLines, 1, 1) == "+")
     Barcodes <- paste("+BX:", substr(ScannedLines[idxRead1 + 1], 1, 16), sep = "")
     NewPlus  <- rep(Barcodes, each = 2)
     ScannedLines[idxPlus] <- NewPlus
  
     # Remove barcodes from read 1 sequences
     ScannedLines[idxRead1 + 1] <- substr(ScannedLines[idxRead1 + 1], 17, 
                                  nchar(ScannedLines[idxRead1 + 1]))
  
     # Get indices for read 1 and read 2 and lines belonging to them
     idx1 <- rep(idxRead1, each = 4) + rep(0:3, length(idxRead1))
     idx2 <- rep(idxRead2, each = 4) + rep(0:3, length(idxRead2))
     writeLines(ScannedLines[idx1], OF1)
     writeLines(ScannedLines[idx2], OF2)
     WriteCounter <- WriteCounter + length(idxRead1)  + length(idxRead2)
     cat("Total of", WriteCounter, " reads written out \n")
   }
}

# Close file connections
close(OF1)
close(OF2)

                   
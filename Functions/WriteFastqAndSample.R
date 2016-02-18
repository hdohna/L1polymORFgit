# General description:
#
#    This function writes a list of reads out as a FastQ file and writes a 
#    sample file to be used by qAlign

# Arguments:
#   
#    ReadList: 
#    Folder: character string indicating the path where the fastq and sample
#       file should be saved
#    FilePrefix: character string indicating the prefix of a file name
#    WriteSample: boolean indicating whether a sample file for qAlign should
#        be written out

# Output:
#   
#    pdf file saved at location specified in Filenamepath

# Comment:

WriteFastqAndSample <- function(ReadList, 
                                Folder,
                                FilePrefix = "Peak",
                                WriteSample = F){
  # Create file name and path for fastq file
  FileNameFastq <- paste(FilePrefix, "fastq", sep = ".")
  FilePathFastq <- paste(Folder, FileNameFastq, sep = "")
  
  # Create lines of the fastq file and save them
  Seqs   <- as.character(ReadList$seq)
  Quals  <- as.character(ReadList$qual)
  SeqIDs <- paste(ReadList$rname, ReadList$pos, sep = ":")
  FastQLines1 <- paste("@", ReadList$qname, sep = "")
  FastQLines3 <- paste("+", SeqIDs, sep = "")
  FastQLinesAll <- rep(NA, 4 * length(FastQLines1))
  FastQLinesAll[seq(1, length(FastQLinesAll) - 3, 4)] <- FastQLines1
  FastQLinesAll[seq(2, length(FastQLinesAll) - 2, 4)] <- Seqs
  FastQLinesAll[seq(3, length(FastQLinesAll) - 1, 4)] <- FastQLines3
  FastQLinesAll[seq(4, length(FastQLinesAll) - 0, 4)] <- Quals
  writeLines(FastQLinesAll, FilePathFastq)
  
  # Create lines of the sample file (for function qAlign) and save them
  FilePathSample <- NULL
  if (WriteSample) {
    FileTable <- rbind(c("FileName",	"SampleName"), 
                       c(FileNameFastq, FilePrefix))
    
    # Create file name and path for fastq file
    FileNameSample <- paste(FilePrefix, "Sample", sep = "")
    FilePathSample <- paste(Folder, FileNameSample, sep = "")
    write.table(FileTable, FilePathSample, sep = "\t", quote = F, 
                row.names = F, col.names = F)
    
  }
  FilePathSample
}




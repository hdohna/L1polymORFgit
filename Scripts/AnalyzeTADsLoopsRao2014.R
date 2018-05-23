# Auxiliary unction to create genomic ranges for both sides of a loop
getLoopGRs <- function(FileName){
  Loops <- read.delim(FileName)
  Loops$chr1 <- paste("chr", Loops$chr1, sep = "")
  Loops$chr2 <- paste("chr", Loops$chr2, sep = "")
  LoopsGR1 <- makeGRangesFromDataFrame(Loops, seqnames.field = "chr1", 
                                       start.field="x1", end.field = "x2")
  LoopsGR2 <- makeGRangesFromDataFrame(Loops, seqnames.field = "chr2", 
                                       start.field="y1", end.field = "y2")
  blnOverlapLoop <- overlapsAny(LoopsGR1, LoopsGR2)
  AllLoops <- c(LoopsGR1, LoopsGR2[!blnOverlapLoop])
  GRangesList(LoopsGR1 = LoopsGR1, LoopsGR2 = LoopsGR2, AllLoops = AllLoops)
}

getLoopDomainGRs <- function(FileName){
  Loops <- read.delim(FileName)
  Loops$chr1 <- paste("chr", Loops$chr1, sep = "")
  Loops$chr2 <- paste("chr", Loops$chr2, sep = "")
  makeGRangesFromDataFrame(Loops, seqnames.field = "chr1", 
                                       start.field="x1", end.field = "y2")
}

# List of files with domains and loops
DomainFiles <- list.files(HiCFolderPath, pattern = "domainlist.txt",
                          full.names = T)
DomainFiles <- DomainFiles[! DomainFiles %in% grep(".gz", DomainFiles, value = T)]
LoopFiles <- list.files(HiCFolderPath, pattern = "looplist.txt",
                        full.names = T)
LoopFiles <- LoopFiles[! LoopFiles %in% grep(".gz", LoopFiles, value = T)]

# Get list of domain and loop ranges
DomainGRList   <- lapply(DomainFiles, getLoopGRs)
LoopGRList     <- lapply(LoopFiles, getLoopGRs)
LoopDomainGRList     <- lapply(LoopFiles, getLoopDomainGRs)
names(DomainGRList) <- sapply(DomainFiles, 
                              function(x) strsplit(x, "_")[[1]][2])
names(LoopGRList) <- sapply(LoopFiles, 
                            function(x) strsplit(x, "_")[[1]][2])
names(LoopDomainGRList) <- sapply(LoopFiles, 
                            function(x) strsplit(x, "_")[[1]][2])

# C
all(names(DomainGRList) == names(LoopGRList))
sapply(1:length(DomainGRList), function(i) {
  sum(overlapsAny(LoopGRList[[i]][[1]], DomainGRList[[i]][[1]]))/length(LoopGRList[[i]][[1]])
})

sapply(1:length(DomainGRList), function(i) {
  sum(overlapsAny(DomainGRList[[i]][[1]], LoopGRList[[i]][[1]]))/length(DomainGRList[[i]][[1]])
})

hist(width(DomainGRList[[1]][[1]]))
hist(width(LoopGRList[[1]][[1]]))

##############################################
#
# General description:
#
#   The following function reads in a table of Histone marks and aggregates
#   them per a set of genomic ranges
#   

# Input:
#
#     FileName: character string providing path to file name with histone marks
#     AggGR: Genomic ranges to aggregate over


# Output:
#   
#     HistAgg: data.frame that contains start and end of AggGR
#          and the mean value of histone

##############################################

AggregateHistoneMarks <- function(FileName, AggGR) {
  
  # Read lines and create table with histone values and genomicRanges
  HistoLines <- readLines(FileName)
  NameSplt   <- strsplit(FileName, "/")[[1]]
  HistName   <- NameSplt[length(NameSplt)]
  HistName   <- strsplit(HistName, "per")[[1]][1]
  blnChr     <- substr(HistoLines, 1, 3) == "chr"
  HistoTable <- read.delim(text = HistoLines[blnChr], header = F,
                           col.names = c("chromosome", "start", "end", HistName))
  HistoGR    <- makeGRangesFromDataFrame(HistoTable)
  
  # Aggregate histone values per range in AggGR   
  Overlaps <- findOverlaps(HistoGR, AggGR)
  HistAgg  <- aggregate(
    HistoTable[Overlaps@from, HistName] ~ Overlaps@to, FUN = mean)
  colnames(HistAgg) <- c("idxGR", HistName)
  HistAgg$start <- start(AggGR)[HistAgg$idxGR]
  HistAgg$end   <- end(AggGR)[HistAgg$idxGR]
  HistAgg$chromosome <- as.vector(seqnames(AggGR))[HistAgg$idxGR]
  HistAgg
}



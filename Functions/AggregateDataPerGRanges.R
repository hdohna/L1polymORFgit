##############################################
#
# General description:
#
#   The following function aggregates a data frame with genomic coordinates 
#   within a set of genomic ranges
#   

# Input:
#
#     RangeWidth:     integer value specifying the width of the range for which
#                     data are summarized
#     SummaryGR:      genomic ranges in which data are summarized. If not 
#                     provided, ranges will be calculated based on RangeWidth 
#                     argument
#     Data2Summarize: data.frame with data to summarize
#     SeqNameCol:     character string of sequence name column name
#     StartCol:       character string of start column
#     EndCol:         character string of end column
#     Cols2Summarize: character vector of column names containing data to 
#                     summarize
#     ChromLengths:   named integer vector of chromosome lengths
#     NoValue: what to enter when no value is provided in a range

# Output:
#   
#     list with follwoing components:
#     SummarizedData: data.frame with summarized data
#     SummaryGR: Genomic ranges for which data were summarized

##############################################

AggregateDataPerGRanges <- function(
  Data2Summarize,
  Cols2Summarize = NULL,
  SeqNameCol = "Chrom",
  StartCol = "genoStart",
  EndCol = "genoEnd",
  RangeWidth = 10^6, 
  SummaryGR = NULL,
  ChromLengths = ChromLengthsHg19,
  SummaryFun = mean,
  NoValue = NA,
  blnAddGRInfo = F,
  blnReturnDFOnly = T
   ){
  
  # Summarize all columns by defaults
  if (is.null(Cols2Summarize)){
    Cols2Summarize <- colnames(Data2Summarize)
  }
  
  # Basic checks that column names are present in input data frame (Data2Summarize)
  if (!SeqNameCol %in% colnames(Data2Summarize)) {
    stop("Column name for sequence names is not input data\n")
  }
  if (!StartCol %in% colnames(Data2Summarize)) {
    stop("Column name for start positions is not input data\n")
  }
  if (!EndCol %in% colnames(Data2Summarize)) {
    stop("Column name for end positions is not input data\n")
  }
  if (any(!Cols2Summarize %in% colnames(Data2Summarize))) {
    stop("Column name for data to summarize is not input data\n")
  }
    
  # Function to create summary genomic ranges
  CreateGR <- function(i, RangeWidth){
    CL <- ChromLengths[i]
    GRStarts <- seq(1, CL, RangeWidth)
    GREnds  <- c(GRStarts[-1] - 1, CL)
    GRanges(seqnames = names(CL), IRanges(start = GRStarts, end = GREnds))
  }
  
  # Create summary genomic ranges if none are supplied
  if (is.null(SummaryGR)){
    # Create genomic ranges for data summaries
    SummaryGR <- CreateGR(1, RangeWidth)
    for (i in 2:length(ChromLengths)){
      NewGR <- CreateGR(i, RangeWidth)
      SummaryGR <- suppressWarnings(c(SummaryGR, NewGR))
    }
  }
  
  # Create genomic ranges for the data to be summarized
  DataGR <- makeGRangesFromDataFrame(Data2Summarize, seqnames.field = SeqNameCol,
                                     start.field = StartCol,
                                     end.field   = EndCol)
  DataSummaryOL <- findOverlaps(DataGR, SummaryGR)

  # Create a data.frame with summarized data
  SummarizedData <- data.frame(matrix(NoValue, nrow = length(SummaryGR),
                                      ncol = length(Cols2Summarize)))
  colnames(SummarizedData) <-  Cols2Summarize
  
  # Calculate mean DNAse per summary genomic range
  SummarizedDataAgg <- aggregate.data.frame(
       Data2Summarize[DataSummaryOL@from, Cols2Summarize], 
       by = list(DataSummaryOL@to), FUN = SummaryFun)
  
  # Replace rows by summaries
  SummarizedData[SummarizedDataAgg$Group.1,] <- 
    SummarizedDataAgg[,seq_along(Cols2Summarize) + 1]
  
  # Add start, end, and chromosome columns
  if(blnAddGRInfo){
    SummarizedData$Start <- start(SummaryGR)
    SummarizedData$End   <- end(SummaryGR)
    SummarizedData$Chrom <- as.vector(seqnames(SummaryGR))
  }
  
  # Return data.frame with summarized data
  if(blnReturnDFOnly) {
    SummarizedData
  
  } else {# Return list with summary genomic ranges
    list(SummarizedData = SummarizedData, SummaryGR = SummaryGR)
  }
}



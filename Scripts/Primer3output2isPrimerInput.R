# The script below reads turns Primer3 output into input for isPrimer

# Read P3 output
P3OutputPath <- "D:/L1polymORF/Data/L1catalog.primer3.output.txt"
P3Output <- readLines(P3OutputPath)

# Get start and end lines for each new element 
NewElementStart <- grep("PRIMER PICKING RESULTS FOR ", P3Output)
NewElementEnd   <- c(NewElementStart[-1] - 1, length(P3Output))

# Define auxilliary function to get primers from primer lines
PrFromLines <- function(PrLines){
  sapply(PrLines, function(x){
    SplitLine <- strsplit(x, " ")[[1]]
    SplitLine[length(SplitLine)]
  })
}

# Define auxilliary function to get product length
LengthFromLines <- function(PrLines){
  sapply(PrLines, function(x){
    SplitLine <- strsplit(x, " ")[[1]]
    SplitLine <- SplitLine[SplitLine != ""]
    LengthNum <- gsub(",", "", SplitLine[3])
    as.numeric(LengthNum)
  })
}

# Define auxilliary function to get info from primer lines
InfoFromLines <- function(PrLines){
  t(sapply(PrLines, function(x){
    SplitLine <- strsplit(x, " ")[[1]]
    SplitLine <- SplitLine[SplitLine != ""]
    StartNum  <- length(grep("[1-9]", SplitLine[1]))
    Info <- SplitLine[c(4:6, 10) + StartNum]
    names(Info) <- c("Length", "Tm", "GC", "Sequence")
    Info
  }))
}

# Loop through elements, extract primers and append them to a growing matrix
OutputLines <- c()
Info <- data.frame()

i <- 1
for (i in 1:length(NewElementStart)){
  idxStart    <- NewElementStart[i]
  idxEnd      <- NewElementEnd[i]
  StartLine   <- P3Output[idxStart]
  LineSubset  <- P3Output[idxStart:idxEnd]
  ElementName <- strsplit(StartLine, "PRIMER PICKING RESULTS FOR ")[[1]][2]
  idxLeftPr   <- grep("LEFT PRIMER", LineSubset)
  idxRightPr  <- grep("RIGHT PRIMER", LineSubset)
  idxLength  <- grep("PRODUCT SIZE:", LineSubset)
  if (length(idxLeftPr) != length(idxRightPr)){
    stop("Not the same number of left and right primers!")
  }
  if (length(idxLeftPr) > 0){
    Names <- paste(ElementName, 1:length(idxLeftPr), sep = "_")
    LeftPrimers  <- PrFromLines(LineSubset[idxLeftPr])
    RightPrimers <- PrFromLines(LineSubset[idxRightPr])
    LeftPrimerInfo  <- InfoFromLines(LineSubset[idxLeftPr])
    RightPrimerInfo <- InfoFromLines(LineSubset[idxRightPr])
    Lengths <- LengthFromLines(LineSubset[idxLength])
    colnames(LeftPrimerInfo) <- paste(colnames(LeftPrimerInfo), "Left", 
                                      sep = "_")
    colnames(RightPrimerInfo) <- paste(colnames(RightPrimerInfo), "Right", 
                                      sep = "_")
    NewPrimerInfo <- cbind(Name = Names, ProductSize = Lengths, LeftPrimerInfo,
                           RightPrimerInfo)
    Info <- rbind(Info, NewPrimerInfo)
    
    #    NewMat       <- cbind(Names, LeftPrimers, RightPrimers)
    NewLines      <- paste(Names, LeftPrimers, RightPrimers, sep = "\t")
    # rownames(NewMat) <- NULL
    # PrimerMat <- rbind(PrimerMat, NewMat)

    OutputLines <- c(OutputLines, NewLines)
  } else {
    warning("No primers for ", ElementName, "\n")
  }
}

#PrimerMat <- PrimerMat[-1, ]
all(OutputLines == paste(Info$Name, Info$Sequence_Left, Info$Sequence_Right, sep = "\t"))

# Write out primer Info
InfoOutputPath <- gsub("primer3.output", "primer3.info", P3OutputPath)[[1]][1]
InfoOutputPath <- gsub(".txt", ".csv", InfoOutputPath)
write.csv(Info, InfoOutputPath, row.names = F)

# Write out left and right primer to be analyzed by isPcr
OutputPath <- gsub("primer3.output", "isPcr.input", P3OutputPath)[[1]][1]
con <- file(OutputPath, "wb")
writeLines(OutputLines, con = con)
close(con)



# Get chromosome length
load("D:/L1polymORF/Data/ChromLengthsHg19.RData")


load("D:/L1polymORF/Data/AnalyzedPacBioL1Ranges_hg19.RData")
ObjectList <- ls()
for (ObjectName in ObjectList){
  NewName <- paste(ObjectName, "PacBio_hg19", sep = "_")
  assign(NewName, eval(parse(text = ObjectName)),
         envir = .GlobalEnv)
}

load("D:/L1polymORF/Data/AnalyzedPacBioL1Ranges_hg19L1.RData")
ObjectList <- ls()
for (ObjectName in ObjectList){
  NewName <- paste(ObjectName, "PacBio_hg19L1", sep = "_")
  assign(NewName, eval(parse(text = ObjectName)),
         envir = .GlobalEnv)
}

load("D:/L1polymORF/Data/NA12878_L1Ranges10X.RData")
for (ObjectName in ObjectList){
  NewName <- paste(ObjectName, "10X", sep = "_")
  assign(NewName, eval(parse(text = ObjectName)),
         envir = .GlobalEnv)
}

max(width(IslGRanges_reduced_10X[idxSuspectL1Ranges_10X]))
OverlapRanges1 <- subsetByOverlaps(IslGRanges_reduced_PacBio_hg19[idxSuspectL1Ranges_PacBio_hg19], 
                                   IslGRanges_reduced_10X[idxSuspectL1Ranges_10X])
OverlapRanges2 <- subsetByOverlaps(IslGRanges_reduced_PacBio_hg19L1[idxSuspectL1Ranges_PacBio_hg19L1], 
                                   IslGRanges_reduced_10X[idxSuspectL1Ranges_10X])
OverlapRanges1B <- subsetByOverlaps(IslGRanges_reduced_10X[idxSuspectL1Ranges_10X],
                                    IslGRanges_reduced_PacBio_hg19[idxSuspectL1Ranges_PacBio_hg19])
length(OverlapRanges1) / length(idxSuspectL1Ranges_PacBio_hg19)
length(OverlapRanges2) / length(idxSuspectL1Ranges_PacBio_hg19L1)

# Get ranges where PacBio reads map
PRanges        <- IslGRanges_reduced_PacBio_hg19[idxSuspectL1Ranges_PacBio_hg19]
PRCountMatched <- table(seqnames(PRanges))
PRCountMatched <- PRCountMatched[names(ChromLengthsHg19)]
PropOverlap <- sapply(1:length(ChromLengthsHg19), function(i){
  ChName   <- names(ChromLengthsHg19)[i]
  Starts   <- sample(ChromLengthsHg19[i], PRCountMatched[i], replace = T)
  PRWidths <- width(PRanges)[as.character(seqnames(PRanges)) == ChName]
  GRSampled <- GRanges(seqnames = ChName, 
                       ranges = IRanges(start = Starts, end = Starts + PRWidths))
  sum(overlapsAny(GRSampled, IslGRanges_reduced_10X[idxSuspectL1Ranges_10X]))/length(GRSampled)
})
mean(PropOverlap)


mean(width(IslGRanges_reduced_PacBio_hg19[idxSuspectL1Ranges_PacBio_hg19]))



# export.bed(OverlapRanges, "D:/L1polymORF/Data/OverlapPacBio10X_hg19")


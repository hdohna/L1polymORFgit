load("D:/L1polymORF/Data/AnalyzedPacBioL1Ranges.RData")
ObjectList <- ls()
for (ObjectName in ObjectList){
  NewName <- paste(ObjectName, "PacBio", sep = "_")
  assign(NewName, eval(parse(text = ObjectName)),
         envir = .GlobalEnv)
}

load("D:/L1polymORF/Data/NA12878_L1Ranges10X.RData")
for (ObjectName in ObjectList){
  NewName <- paste(ObjectName, "10X", sep = "_")
  assign(NewName, eval(parse(text = ObjectName)),
         envir = .GlobalEnv)
}
IslGRanges_reduced_10X
IslGRanges_reduced_PacBio

OverlapRanges <- subsetByOverlaps(IslGRanges_reduced_PacBio[idxSuspectL1Ranges_PacBio], 
                 IslGRanges_reduced_10X[idxSuspectL1Ranges_10X])
length(idxSuspectL1Ranges_PacBio)
length(idxSuspectL1Ranges_10X)

# export.bed(OverlapRanges, "D:/L1polymORF/Data/OverlapPacBio10X_hg19")


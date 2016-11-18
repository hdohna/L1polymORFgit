# The following script compares results from PacBio and PCR

# read result from pacBio (produced by scripts AlignLocalL1Peaks and 
# AnalyzeL1PCResults)
PacBioResults <- read.csv("D:/L1polymORF/Data/L1PacBioResults.csv")
PCRResults    <- read.csv("D:/L1polymORF/Data/NA12878_L1PCR_results_filtered.csv")
colnames(PCRResults)[colnames(PCRResults) == "AccNr"] <- "Accession"

# Merge results
MergedResults <- merge(PacBioResults, PCRResults, by = "Accession")
L1PresentPacBio <- MergedResults$MaxCover > 0
L1PresentPCR    <- MergedResults$ProbPresent > 0.95

# Create a table that allows comparing results
ResultComparison <- MergedResults[,c("Accession", "ProbPresent", "MaxCover", 
                                     "chrom", "start")]

# Calculate contingency table to 
PCR_PacBioComparison <- table(L1PresentPacBio, L1PresentPCR)
dimnames(PCR_PacBioComparison)
write.table(PCR_PacBioComparison, "D:/L1polymORF/Data/PCR_PacBioComparison")
chisq.test(L1PresentPacBio, L1PresentPCR)
sum(L1PresentPacBio & L1PresentPCR & (nchar(as.character(MergedResults$chrom)) > 5))
sum(L1PresentPacBio & L1PresentPCR & (nchar(as.character(MergedResults$chrom)) <= 5))


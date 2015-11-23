# The following script aligns reads from a chip-seq experiment to detect L1 insertions

# Load package
library(QuasR)

# Run alignment for hg19
qAlign(sampleFile = "D:/L1polymORF/Data/SampleFile_NA12878_NA12892.txt",
  genome = "BSgenome.Hsapiens.UCSC.hg19", paired = "fr")

# Run alignment for hg38
AlignHg38 <- qAlign(sampleFile = "D:/L1polymORF/Data/SampleFile_NA12878_NA12892_hg38.txt",
       genome = "BSgenome.Hsapiens.UCSC.hg38", paired = "fr")

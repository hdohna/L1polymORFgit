# The script below runs the function 'driverL1Analysis'

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Load packages
library(ShortRead)
library(csaw)

# Run function
driverL1Analysis(
  PeakBam = "/home/hzudohna/NA12878PacBio_alnMappedSorted.bam", 
  L1HSBamFile = "/home/hzudohna/L1polymORF/Data/NA12878PacBio_L1hg19.bam", 
  FastQFolder = NULL, 
  L1HSConsensus = "/home/hzudohna/L1polymORF/Data/Homo_sapiens_L1_consensus.fa",
  L1RefRanges    = '/home/hzudohna/L1polymORF/Data/L1RefRanges_hg19.Rdata',
  OutputFolder = "/home/hzudohna/L1polymORF/Data/", 
  CoverSummaryPlot, 
  CoverComparePlot, 
  ResultFileName = "L1HSPacBio_Results.Rdata",
  MinMaxCover = 5, 
  MinGap = 5, 
  MinDist2L1 = 20000, 
  blnComparePeaksWithRefL1 = T,
  blnWriteFastq     = F,
  blnMap2L1         = F, 
  blnAddReadGroups  = T, 
  blnCreateBamIndices = F,
  blnCallHaplotypes = T, 
  blnAnalyze        = F,
  AlignCommand = '/home/txw/bwa/bwa-0.7.12/bwa mem',
  ReadGroupSuffix = "withRG.bam",
  BamSuffix = "withRG.bam",
  HapTypeCallOptions = "",
  BamSuffixHapTypeCall = "withRG.bam"
)

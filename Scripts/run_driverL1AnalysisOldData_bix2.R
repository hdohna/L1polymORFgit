# The script below runs the function 'driverL1Analysis'

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Load packages
library(ShortRead)
library(csaw)

# Run function
driverL1Analysis(
  PeakBam = "/home/hzudohna/BoData/NA12878-L1HS_S1_L001.dedup.unique.sorted.bam", 
  L1HSBamFile = NULL, 
  FastQFolder = "/home/hzudohna/L1polymORF/Data/FastQ", 
  L1HSConsensus = "/home/hzudohna/L1polymORF/Data/Homo_sapiens_L1_consensus.fa",
  L1RefRanges    = '/home/hzudohna/L1polymORF/Data/L1RefRanges_hg38.Rdata',
  OutputFolder = "/home/hzudohna/L1polymORF/Data/", 
  CoverSummaryPlot, 
  CoverComparePlot, 
  ResultFileName = "NA12878-L1HS_Results.Rdata",
  MinMaxCover = 5, 
  MinGap = 6000, 
  MinDist2L1 = 20000, 
  blnComparePeaksWithRefL1 = T,
  blnWriteFastq     = T,
  blnMap2L1         = T, 
  blnAddReadGroups  = T, 
  blnCreateBamIndices = T,
  blnCallHaplotypes = T, 
  blnAnalyze        = F,
  AlignCommand = '/home/txw/bwa/bwa-0.7.12/bwa mem',
  ReadGroupSuffix = "withRG.bam",
  BamSuffix = "withRG.bam",
  HapTypeCallOptions = "",
  BamSuffixHapTypeCall = "withRG.bam"
)

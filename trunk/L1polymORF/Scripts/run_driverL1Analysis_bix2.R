# The script below runs the function 'driverL1Analysis'

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Load packages
library(ShortRead)
library(csaw)

# Run function
driverL1Analysis(
  PeakBam = "/home/hzudohna/NA12878-L15P_S1_L001_001.dedup.mapnonzero.sorted.bam", 
#  PeakBam = "/share/diskarray2/L1HS/fastq_and_align/LHS0001comb_hg19_sorted_rmdupMKDP_rln_recal.bam", 
  L1HSBamFile = NULL, 
  FastQFolder = "/share/diskarray2/L1HS/fastq_and_align/LHS0001combtrimRCNMP3NMP3Q20PE-33564531/", 
  L1HSConsensus = "/home/hzudohna/L1polymORF/Data/Homo_sapiens_L1_consensus.fa",
  L1RefRanges    = '/home/hzudohna/L1polymORF/Data/L1RefRanges_hg19.Rdata',
  OutputFolder = "/home/hzudohna/L1polymORF/Data/", 
  CoverSummaryPlot, 
  CoverComparePlot, 
ResultFileName = "LHS0001combHG38_Results.Rdata",
#ResultFileName = "LHS0001comb_Results.Rdata",
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

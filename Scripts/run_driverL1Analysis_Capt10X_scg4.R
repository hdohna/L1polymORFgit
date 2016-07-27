# The script below runs the function 'driverL1Analysis'

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load packages
library(ShortRead)
library(csaw)

# Run function
driverL1Analysis(
  PeakBam = "/srv/gsfs0/projects/levinson/hzudohna/10Xcapture/NA12878_capt10X_aln2hg19.bam", 
  L1HSBamFile = NULL, 
  FastQFolder = "/srv/gsfs0/projects/levinson/hzudohna/10Xcapture/FastqFiles/", 
  L1HSConsensus = "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/Homo_sapiens_L1_consensus.fas",
  L1RefRanges    = '/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/L1RefRanges_hg19.Rdata',
  OutputFolder = "/srv/gsfs0/projects/levinson/hzudohna/10Xcapture/", 
  ResultFileName = "NA12878_capt10X_Results.Rdata",
  PlotFolder = "/srv/gsfs0/projects/levinson/hzudohna/10Xcapture/Figures/", 
  MinMaxCover = 1, 
  MinGap = 1000, 
  MinDist2L1 = 20000, 
  NrChromPieces = 20,
  blnComparePeaksWithRefL1 = T,
  blnWriteFastq     = T,
  blnMap2L1         = T, 
  blnAddReadGroups  = T, 
  blnCreateBamIndices = T,
  blnCallHaplotypes = F, 
  blnAnalyze        =   T,
  IdChar2Remove = 4,
  AlignCommand = c('module load bwa', 'bwa mem'),
  IndexCommand = c('module load bwa', 'bwa index'),
  ReadGroupSuffix = "withRG.bam",
  BamSuffix = "withRG.bam",
  HapTypeCallOptions = "",
  BamSuffixHapTypeCall = "withRG.bam"
)

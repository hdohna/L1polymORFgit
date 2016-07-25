# The script below runs the function 'driverL1Analysis'

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load packages
library(ShortRead)
library(csaw)

# Run function
driverL1Analysis(
  PeakBam = "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/BZ_NA12878L1capt5-9kb_readsofinsert_hg19.bam", 
  L1HSBamFile = NULL, 
  FastQFolder = "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/FastqFiles/", 
  L1HSConsensus = "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/Homo_sapiens_L1_consensus.fas",
  L1RefRanges    = '/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/L1RefRanges_hg19.Rdata',
  OutputFolder = "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/", 
  ResultFileName = "BZ_NA12878L1capt5-9kb_Results.Rdata",
  PlotFolder = "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/Figures/", 
  MinMaxCover = 1, 
  MinGap = 6000, 
  MinDist2L1 = 20000, 
  blnComparePeaksWithRefL1 = F,
  blnWriteFastq     = F,
  blnMap2L1         = F, 
  blnAddReadGroups  = F, 
  blnCreateBamIndices = T,
  blnCallHaplotypes = F, 
  blnAnalyze        =   T,
  IdChar2Remove = 4,
  AlignCommand = c('module load bwa', 'bwa mem -k17 -W40 -r10 -A2 -B5 -O2 -E1 -L0'),
  IndexCommand = c('module load bwa', 'bwa index'),
  ReadGroupSuffix = "withRG.bam",
  BamSuffix = "withRG.bam",
  HapTypeCallOptions = "",
  BamSuffixHapTypeCall = "withRG.bam"
)

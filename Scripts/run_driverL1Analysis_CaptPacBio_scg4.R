# The script below runs the function 'driverL1Analysis'

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load packages
library(ShortRead)
library(csaw)

# Run function
driverL1Analysis(
  PeakBam = "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/PacBioL1captPS/PacBioHiFiSubreads_hg19masked_ngmlr.sorted.bam", 
#  L1HSBamFile = "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/BZ_NA12878L1capt5-9kb_sub_reads_aln2L1.sorted.bam", 
  L1HSConsensus = "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/Homo_sapiens_L1_consensus.fas",
  L1RefRanges    = '/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/L1RefRanges_hg19.Rdata',
  OutputFolder = "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/", 
  ResultFileName = "PacBioHiFi__Results.Rdata",
  PlotFolder = "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/Figures/", 
  MinMaxCover = 2, 
  MinGap = 100, 
  MinDist2L1 = 500, 
  NrChromPieces = 1,
  blnComparePeaksWithRefL1 = T,
  blnFilterOverlap         = F,
  blnWriteFastq            = T,
  blnMap2L1                = T, 
  blnSam2Bam               = T,
  blnCreateBamIndices      = T,
  blnCalcCoverMat          = T,
  IdChar2Remove = 4,
  NrJobsPerBatch = 50, 
  WaitBetwJobs = 60,
  NrReadsPerIter = 10^3,
  AlignCommand = c('module load bwa', 'bwa mem -k17 -W40 -r10 -A2 -B5 -O2 -E1 -L0'),
  IndexCommand = c('module load bwa', 'bwa index'),
  BamSuffix = ".bam"
)

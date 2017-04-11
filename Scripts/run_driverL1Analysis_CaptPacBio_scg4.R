# The script below runs the function 'driverL1Analysis'

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load packages
library(ShortRead)
library(csaw)

# Run function
driverL1Analysis(
  PeakBam = "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/BZ_NA12878L1capt5-9kb_subreads_hg19masked.sorted.bam", 
#  L1HSBamFile = "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/BZ_NA12878L1capt5-9kb_sub_reads_aln2L1.sorted.bam", 
  FastQFolder = "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/FastqSubreads/", 
  L1HSConsensus = "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/Homo_sapiens_L1_consensus.fas",
  L1RefRanges    = '/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/L1RefRanges_hg19.Rdata',
  OutputFolder = "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/", 
  ResultFileName = "BZ_NA12878L1capt5-9kb_Results.Rdata",
  PlotFolder = "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/Figures/", 
  MinMaxCover = 1, 
  MinGap = 100, 
  MinDist2L1 = 2000, 
  NrChromPieces = 1,
  blnComparePeaksWithRefL1 = F,
  blnWriteFastq            = F,
  blnBam2Fastq             = F,
  blnMap2L1                = F, 
  blnSam2Bam        = F,
  blnCreateBamIndices      = F,
  blnCalcCoverMat          = T,
  IdChar2Remove = 4,
  NrJobsPerBatch = 50, 
  WaitBetwJobs = 100,
  NrReadsPerIter = 10^3,
  AlignCommand = c('module load bwa', 'bwa mem -k17 -W40 -r10 -A2 -B5 -O2 -E1 -L0'),
  IndexCommand = c('module load bwa', 'bwa index'),
  BamSuffix = ".bam"
)

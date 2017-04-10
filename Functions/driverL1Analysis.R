##############################################
#
# General description:
#
#   The following function reads in a bam file that contains peaks around L1HS
#   insertion sites. It creates per peak a bam file of reads mapped to L1HS and
#   a vcf file.

# Input:
#
#  Folders and file names:
#     PeakBam:       bam file containing peaks around L1HS
#     L1HSBamFile:   bam file with reads mapped to L1HS
#     FastQFolder:   Folder containing the original fast
#     L1HSConsensus: fasta file containing the L1HS consensus sequence
#     OutputFolder:  folder where analysis data output is saved
#     CoverSummaryPlot: path where pdf file with plot of coverage summary is 
#                       saved
#     CoverComparePlot:path where pdf file with plot of coverage summary is 
#                       saved

#  Peak calling parameters:
#     MinMaxCover: minimum maximum coverage to be called a peak 
#     MinGap: minimum gap allowed between two separate peaks. Peaks with 
#        smaller gap than this will be collapsed into one
#     MinDist2L1: minimum distance to reference L1 to be considered 
#        "non-reference" 
#     PacBioWindow <- 1000

#  Run parameters:
#     blnGetL1Ranges: boolean indicator whether L1 ranges should be 
#        determined
#     blnWriteFastq: boolean indicator whether fastq files per peak should be
#     blnMap2L1: boolean indicator whether 
#     blnAddReadGroups: boolean indicator whether 
#     blnCallHaplotypes: boolean indicator whether 
#     blnAnalyze: boolean indicator whether 

#  Run commands:
#     AlignCommand: character string with command to run alignment to reference
#         L1HS
#     AddGroupCmd: character string with command to add read groups to bam 
#         files
#     AddGroupOptions: character vector specifying options add read groups to
#         bam files

# Output:
#   
#    Little sam file for each fastq file

##############################################

driverL1Analysis <- function(
  PeakBam, 
  L1HSBamFile = NULL, 
  FastQFolder, 
  L1HSConsensus = "/home/hzudohna/L1polymORF/Data/Homo_sapiens_L1_consensus.fa",
  L1RefRanges    = '/home/hzudohna/L1polymORF/Data/L1RefRanges_hg19.Rdata',
  OutputFolder = "/home/hzudohna/L1polymORF/Data/", 
  PlotFolder = "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/Figures/", 
  ResultFileName,
  MinMaxCover, MinGap, 
  MinDist2L1, 
  NrChromPieces = 1,
  blnComparePeaksWithRefL1 = F,
  blnWriteFastq     = F,
  blnFilterBamPerL1 = F,
  blnBam2Fastq = F,
  blnMap2L1         = F, 
  blnAddReadGroups  = F,
  blnCreateBamIndices = F,
  blnFilterBam = F,
  blnCallHaplotypes = F, 
  blnCalcCoverMat = F,
  IdChar2Remove = 4,
  EndList = NULL,
  blnFilterOverlap = F,
  NrJobsPerBatch = 100, 
  WaitBetwJobs = 1000,
  NrReadsPerIter = 10^6,
  AlignCommand = c('module load bwa', 'bwa mem'),
  IndexCommand = c('module load bwa', 'bwa index'),
  AddGroupCmd  = c('module load picard-tools/2.0.1', 
                   "java -jar picard.jar AddOrReplaceReadGroups"),
  AddGroupOptions = c("RGLB=lib1", "RGPL=illumina", "RGPU=unit1", "RGSM=20",
                      "SORT_ORDER=null", "CREATE_INDEX=TRUE", "VALIDATION_STRINGENCY=LENIENT"),
  CreateBamIndexCmd = c('module load samtools', 'samtools index'),
  HapTypeCallCmd = c("module load gatk/3.4.0", 
                     "java -jar GenomeAnalysisTK.jar -T HaplotypeCaller"),
  HapTypeCallOptions = "--emitRefConfidence GVCF",
  BamSuffix = ".bam",
  ReadGroupSuffix = "withRG.bam",
  BamSuffixHapTypeCall = ".bam",
  SamSuffix = ".sam",
  VCFSuffix = ".vcf" 
){
  
  #######################################################
  #                                                     #
  #          Create folders and file paths              #
  #                                                     #
  #######################################################
  
  # Specify file paths 
  InFileSplit <- strsplit(PeakBam, "/")[[1]]
  InFileName  <- InFileSplit[length(InFileSplit)]
  FolderPrefix <- strsplit(InFileName, "_")[[1]][1]
  
  # Create new folder names and folders
  OutFolderName_NonRef <- paste(OutputFolder, FolderPrefix, "_NonRef", 
                                sep = "")
  OutFolderName_Ref    <- paste(OutputFolder, FolderPrefix, "_Ref", sep = "")
  if (!dir.exists(OutFolderName_NonRef)) {
    dir.create(OutFolderName_NonRef)
    cat("Created folder", OutFolderName_NonRef, "\n")
  }
  if (!dir.exists(OutFolderName_Ref)) {
    dir.create(OutFolderName_Ref)
    cat("Created folder", OutFolderName_Ref, "\n")
  }
  
  # Create path to result file
  OutResults_RangeComparison  <- paste(OutputFolder, FolderPrefix, 
                                       "_L1Ranges.RData", sep = "")
  OutResults_Analysis  <- paste(OutputFolder, ResultFileName, sep = "")
  
  # Create paths for plots
  CoverSummaryPlot <- paste(PlotFolder, FolderPrefix, 
                            "_CoverSummary.pdf", sep = "")
  CoverComparePlot <- paste(PlotFolder, FolderPrefix, 
                            "_CoverCompare.pdf", sep = "")
  
  #######################################################
  #                                                     #
  #    Compare peaks and reference L1 ranges            #
  #                                                     #
  #######################################################
  
  if(blnComparePeaksWithRefL1){
    
    # Create a name for bam file for reads that map to full-length L1 on the
    # reference genome
    L1RefBamFile <- paste(FolderPrefix, "_L1Ref.bam", sep = "")
    OutBamFileFullLengthL1 <- paste(OutFolderName_Ref, L1RefBamFile, sep = "/")
    ComparePeaksWithRefL1(
      BamFile = PeakBam,
      OutBamFileFullLengthL1 = OutBamFileFullLengthL1,
      L1Ranges = L1RefRanges,
      MinMaxCover = MinMaxCover,    # minimum maximum coverage to be called a peak 
      MinGap      = MinGap,
      MinDist2L1  = MinDist2L1, # minimum distance to L1 to be called a peak 
      NrChromPieces = NrChromPieces,
      OutFile = OutResults_RangeComparison,
      EndList = EndList,
      blnFilterOverlap = blnFilterOverlap)
  }
  
  # Load the results produced by the function ComparePeaksWithRefL1. The main
  # objects contained in this file that are used by subsequent analysis steps
  # are SuspectL1Ranges and idxSuspectL1Ranges. Check whether they are there
  load(OutResults_RangeComparison)
  if (!exists("SuspectL1Ranges")){
    stop("Object SuspectL1Ranges missing!")
  }
  if (!exists("idxSuspectL1Ranges")){
    stop("Object idxSuspectL1Ranges missing!")
  }
  
  # Create names for fastq and bam files that are saved per suspected 
  # non-reference L1
  FilePrefixes     <- paste(seqnames(SuspectL1Ranges), idxSuspectL1Ranges, 
                            sep = "_")
  LittleFastqFiles <- paste(FilePrefixes, ".fastq", sep = "")
  LittleBamFiles   <- paste(FilePrefixes, ".bam", sep = "")
  LittleFastqPaths <- paste(OutFolderName_NonRef, LittleFastqFiles, sep = "/")
  LittleBamPaths   <- paste(OutFolderName_NonRef, LittleBamFiles, sep = "/")
  
  #######################################################
  #                                                     #
  #    Write fastq of suspected L1 not in reference     #
  #                                                     #
  #######################################################
  
  if(blnWriteFastq & is.null(L1HSBamFile)){
    
    # Write little fastq files per suspected peak
    WriteFastQPerRange(Ranges = SuspectL1Ranges, 
                       InBamfilePath  = PeakBam,
                       InFastQfilePaths = list.files(FastQFolder, full.names = T),
                       OutFilePaths = LittleFastqPaths,
                       NrReadsPerIter = NrReadsPerIter,
                       IdChar2Remove = IdChar2Remove) 
  }
  
  OutBamFilePaths <- gsub(".fastq", ".bam", LittleFastqPaths)
  if(blnFilterBamPerL1 & is.null(L1HSBamFile)){
    InBamFilePaths <- FilterBamPerRange(Ranges = SuspectL1Ranges, 
                          InBamfilePath = PeakBam, 
                          OutBamFilePaths = OutBamFilePaths) 
  } 
  if(blnBam2Fastq & is.null(L1HSBamFile)){
    ConvertMultiBam2Fastq(OutBamFilePaths, 
                          NrJobsPerBatch = NrJobsPerBatch, 
                          WaitBetwJobs = WaitBetwJobs) 
  }
  
  #######################################
  #                                     #
  #     Map fastq file per range        #
  #            to L1HS                  #
  #                                     #
  #######################################
  
  if(blnMap2L1 & is.null(L1HSBamFile)){
    
    FilePaths <- MapMultiFastq(FastQFolder  = OutFolderName_NonRef,
                               AlignCommand = AlignCommand,
                               IndexCommand = IndexCommand,
                               Reference = L1HSConsensus,
                               SamSuffix = SamSuffix,
                               NrJobsPerBatch = NrJobsPerBatch, 
                               WaitBetwJobs = WaitBetwJobs)
  }
  
  #######################################
  #                                     #
  #     Turn sam files into bam files   #
  #                                     #
  #######################################
  
  if(is.null(L1HSBamFile)){
    
    # Get all names of sam files created by BWA
    SamFileNames <- list.files(OutFolderName_NonRef, pattern = SamSuffix,
                               full.names = T)
    
    # Turn sam files into bam files
    for (fn in SamFileNames) {
      cat("Turning", fn, "into a bam file\n")
      asBam(fn, destination = substr(fn, 1, nchar(fn) - 4), overwrite = T)
    }
    
  }
  
  ###################################################
  #                                                 #
  #     Filter bam file with reads mapped to L1HS   #
  #                                                 #
  ###################################################
  
  if(!is.null(L1HSBamFile) & blnFilterBam){
    
    # Loop over peaks that do not overlap with reference L1 and write out a
    # per peak a separate bam file of reads mapped to L1HS
    for (i in 1:length(SuspectL1Ranges)){
      param <- ScanBamParam(which = SuspectL1Ranges[i], what = "qname")
      IDs   <- scanBam(PeakBam, param = param)
      IDs   <- unlist(IDs)
      IDFilter <- FilterRules(getIDs <- function(DF){DF$qname %in% IDs})
      cat("Writing filtered L1Hs bam file", LittleBamPaths[i], "\n")
      filterBam(L1HSBamFile, LittleBamPaths[i], filter = IDFilter)
    }
  }
  
  #######################################
  #                                     #
  #     Add read groups                 #
  #                                     #
  #######################################
  
  if(blnAddReadGroups){
    
    BamSuffix   <- ReadGroupSuffix
    BamSuffixHapTypeCall <- ReadGroupSuffix
    FilePathsRG <- AddMultiReadGroups(BamFolder = OutFolderName_NonRef,
                                      AddGroupCmd   = AddGroupCmd,
                                      AddGroupOptions = AddGroupOptions,
                                      ReadGroupSuffix = ReadGroupSuffix)
  }
  
  #######################################
  #                                     #
  #     Create bam indices              #
  #                                     #
  #######################################
  
  if(blnCreateBamIndices){
    
    CreateMultiBamIndex(BamFolder  = OutFolderName_NonRef, CreateBamIndexCmd = CreateBamIndexCmd, 
                        BamSuffix = BamSuffix)
  }
  
  
  
  if(blnCallHaplotypes){
    
    FilePathsVCF <- CallMultiVariants(BamFolder = OutFolderName_NonRef,  
                                      HapTypeCallCmd = HapTypeCallCmd,
                                      RefSeqPath  = L1HSConsensus,
                                      OptionLines = HapTypeCallOptions,
                                      BamSuffix = BamSuffixHapTypeCall,
                                      VCFSuffix = VCFSuffix) 
  }
  
  #######################################
  #                                     #
  #     Calculate coverage on L1        #
  #                                     #
  #######################################
  
  if (blnCalcCoverMat){
    CalcCoverMatReadList(
      OutputFilePath = OutResults_Analysis,
      OutFolderName_NonRef = OutFolderName_NonRef,
      L1RangesPath = OutResults_RangeComparison)
  }
}
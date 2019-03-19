# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load necessary packages
library(Rsamtools)
library(rtracklayer)

# Run parameters
RunTime_MELT <- '12:00:00'
Mem_MELT     <- '500G'
JobsTotal    <- 50

# Path to reference data and for 1000 genomes
#RefPath   <- "D:/L1polymORF/Data/"
RefPath     <- "/labs/dflev/hzudohna/RefSeqData/"
RefFilePath <- "/labs/dflev/hzudohna/RefSeqData/hg19.fa"
RefFilePath <- "/labs/dflev/hzudohna/RefSeqData/hs37d5.fa.gz"

Path1000G   <-   "/labs/dflev/hzudohna/1000Genomes/"

# Path for scratch storage
ScratchPath <- "/scratch/users/hzudohna/"

# Load necessary objects
load(paste(Path1000G, 'GRanges_L1_1000Genomes.RData', sep = ""))

# Create bed file of ranges around each 1000 genome L1
# ChrNames <- as.vector(seqnames(L1_1000G_GR_hg19))
# ChrNrs   <- substr(ChrNames, 4, nchar(ChrNames))
# L1NeighborRanges_1000G <- GRanges(seqnames = ChrNrs, 
#                                   IRanges(start = start(L1_1000G_GR_hg19) - 500,
#                                           end   = end(L1_1000G_GR_hg19) + 500))
# L1NeighborbedPath_1000G <- paste(RefPath, "L1NeighborRanges_1000G.bed", sep = "")
# L1bedPath_1000G         <- paste(RefPath, "L1Ranges_1000G.bed", sep = "")
# export.bed(c(L1Ranges_shortName, L1NeighborRanges_1000G), L1NeighborbedPath_1000G) 

# Specify the general path to 1000 genome bam file
BamPath1000G_General <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/IndividualID/alignment/IndividualID.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam"
DirPath1000G_General <- "ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/IndividualID/alignment/"

# Get all directories in the 1000 genome path
# Dirs1000G <- list.dirs(Path1000G, recursive = F, full.names = F)
# blnNotAnalyzed <- !SampleColumns %in% Dirs1000G
# IndividualID <- SampleColumns[blnNotAnalyzed][3]
# sum(blnNotAnalyzed)

# Get names of vcf files
VcfFiles <- paste("/labs/dflev/hzudohna/1000Genomes/", SampleColumns, "_fullGenome",
                  "/LINE1.final_comp.vcf", sep = "")
blnNoResults <- !file.exists(VcfFiles)
cat(sum(!blnNoResults), "genomes with results\n")
cat(sum(blnNoResults), "genomes without results\n")

# Remove the files in folders with nor results
NoResultFolders <- gsub("/LINE1.final_comp.vcf", "", VcfFiles[blnNoResults])

RemoveCommand <- paste("rm -r", paste(NoResultFolders, collapse = " "))
cat("Removing", length(NoResultFolders))
system(RemoveCommand)
blnNotAnalyzed <- !dir.exists(paste("/labs/dflev/hzudohna/1000Genomes/", 
                                    SampleColumns, "_fullGenome", sep = ""))
cat(sum(!blnNotAnalyzed), "genomes analyzed\n")
cat(sum(blnNotAnalyzed), "genomes not analyzed\n")
SampleColumns[blnNoResults][1:5]

########################################################################################
# Run MELT for low coverage bam files
IndividualID = SampleColumns[blnNotAnalyzed]
IDs2Analyze <- SampleColumns[blnNoResults][1:min(sum(blnNoResults), JobsTotal)]
for (IndividualID in IDs2Analyze){

  # Define paths
  DirPath1000G <- gsub("IndividualID", IndividualID, DirPath1000G_General)
  BamOutPath   <- paste(Path1000G, "L1Filtered_", IndividualID, ".bam", sep = "")
  BamOutPath_Lftp.bam   <- paste(ScratchPath, IndividualID, ".bam", sep = "")
  BamOutPath_Lftp.bai   <- paste(ScratchPath, IndividualID, ".bai", sep = "")
  
  # Create folder for MELT results
  NewDir <- paste(Path1000G, IndividualID, "_fullGenome", sep = "")
  MkDirCmd <- NULL
  if (!dir.exists(NewDir)){
    paste("mkdir", NewDir)
  }
  
  # Create curl command to get file names on ftp directories
  DirListPath  <- paste(Path1000G, IndividualID, "_FileList.txt", sep = "")
  CurlCmd <- paste("curl", paste("ftp://", DirPath1000G, sep = ""),
                   ">", DirListPath)
  system(CurlCmd)

  # Get name of bam file 
  DirListLines <- readLines(DirListPath)
  FileNames <- sapply(DirListLines, function(x){
    SplitLines <- strsplit(x, " ")[[1]]
    SplitLines[length(SplitLines)]
  })
  FileNames   <- FileNames[grep("\\.mapped", FileNames)]
  BamFileName <- FileNames[-grep("bam.", FileNames)]
  BaiFileName <- FileNames[grep(".bai", FileNames)]
  BamPath1000G <- paste("ftp://", DirPath1000G, BamFileName, sep = "")
  BaiPath1000G <- paste("ftp://", DirPath1000G, BaiFileName, sep = "")
  
  # # Construct samtools command to get filtered bam file
  # SamToolsCmds <- c("module load samtools",
  #                   paste("samtools view -b -h", BamPath1000G, "-L", 
  #                         L1NeighborbedPath_1000G, "| samtools sort -o",
  #                         BamOutPath),
  #                   paste("samtools index", BamOutPath))
  
  # Construct lftp command to get  bam file
  LftpCmds <- paste("lftp -c'", 
                    paste("open ftp://", DirPath1000G, ";", sep = ""),
                    paste("get1 -o", BamOutPath_Lftp.bam, " ", BamPath1000G, ";", sep = ""),
                    paste("get1 -o", BamOutPath_Lftp.bai, " ", BaiPath1000G, ";", sep = ""),
                    "'")
  
  # Construct MELT command
  MELTCmds <- c(MkDirCmd,
                "module load bowtie2",
                paste("java -Xmx2g -jar /labs/dflev/hzudohna/MELTv2.1.5/MELT.jar Single -bamfile",
                      BamOutPath_Lftp.bam, 
                      "-c 7",
                      "-h", RefFilePath,
                      "-t /labs/dflev/hzudohna/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip",
                      "-n /labs/dflev/hzudohna/MELTv2.1.5/add_bed_files/1KGP_Hg19/hg19.genes.bed",
                      "-b MT/GL000207.1/GL000226.1/GL000229.1/GL000231.1/GL000210.1/GL000239.1/GL000235.1/GL000201.1/GL000247.1/GL000245.1/GL000197.1/GL000203.1/GL000246.1/GL000249.1/GL000196.1/GL000248.1/GL000244.1/GL000238.1/GL000202.1/GL000234.1/GL000232.1/GL000206.1/GL000240.1/GL000236.1/GL000241.1/GL000243.1/GL000242.1/GL000230.1/GL000237.1/GL000233.1/GL000204.1/GL000198.1/GL000208.1/GL000191.1/GL000227.1/GL000228.1/GL000214.1/GL000221.1/GL000209.1/GL000218.1/GL000220.1/GL000213.1/GL000211.1/GL000199.1/GL000217.1/GL000216.1/GL000215.1/GL000205.1/GL000219.1/GL000224.1/GL000223.1/GL000195.1/GL000212.1/GL000222.1/GL000200.1/GL000193.1/GL000194.1/GL000225.1/GL000192.1/NC_007605",
                      "-w", NewDir))

  # Command to remove bam files
#  RemCmds <- c(paste("rm", BamOutPath), paste("rm ", BamOutPath, ".bai", sep = ""))
  RemCmds <- c(paste("rm", BamOutPath_Lftp.bam), paste("rm ", BamOutPath_Lftp.bai, sep = ""))
  
  # Create script name and launch script
  ScriptName <- paste("ME_L1Ins_Script", IndividualID, sep = "_")
  CreateAndCallSlurmScript(file = ScriptName, 
                           SlurmCommandLines = c(LftpCmds, 
                                                 "echo 'bam file retrieved'", 
                                                 "module load bowtie",
                                                 MELTCmds, 
                                                 "echo 'MELT completed'",
                                                 RemCmds,
                                                 "echo 'bam files removed'"),
                           RunTime = RunTime_MELT,
                           Mem = Mem_MELT)

  
}

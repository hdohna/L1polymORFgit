# The following script creates objects and commands for calling L1 deletions in
# 1000 genome data using MELT.

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load necessary packages
library(Rsamtools)
library(rtracklayer)

# Path to reference data and for 1000 genomes
#RefPath   <- "D:/L1polymORF/Data/"
RefPath   <- "/labs/dflev/hzudohna/RefSeqData/"
Path1000G <- "/labs/dflev/hzudohna/1000Genomes/"

# Load necessary objects
load(paste(Path1000G, 'GRanges_L1_1000Genomes.RData', sep = ""))
load(paste(RefPath, 'L1RefRanges_hg19.Rdata', sep = ""))

# Create bed file of ranges around each reference L1
ChrNames <- as.vector(seqnames(L1GRanges))
ChrNrs   <- substr(ChrNames, 4, nchar(ChrNames))
L1NeighborRanges <- GRanges(seqnames = ChrNrs, 
                            IRanges(start = start(L1GRanges) - 500,
                                    end   = end(L1GRanges) + 500))
L1NeighborbedPath <- paste(RefPath, "L1NeighborRanges.bed", sep = "")
L1bedPath         <- paste(RefPath, "L1Ranges.bed", sep = "")
export.bed(L1NeighborRanges, L1NeighborbedPath) 

L1Ranges_shortName <- GRanges(seqnames = ChrNrs, 
                              IRanges(start = start(L1GRanges),
                                      end   = end(L1GRanges)))
export.bed(L1Ranges_shortName, L1bedPath) 


# Create bed file of ranges around each 1000 genome L1
ChrNames <- as.vector(seqnames(L1_1000G_GR_hg19))
ChrNrs   <- substr(ChrNames, 4, nchar(ChrNames))
L1NeighborRanges_1000G <- GRanges(seqnames = ChrNrs, 
                            IRanges(start = start(L1_1000G_GR_hg19) - 500,
                                    end   = end(L1_1000G_GR_hg19) + 500))
L1NeighborbedPath_1000G <- paste(RefPath, "L1NeighborRanges_1000G.bed", sep = "")
L1bedPath_1000G         <- paste(RefPath, "L1Ranges_1000G.bed", sep = "")
#export.bed(L1NeighborRanges_1000G, L1NeighborbedPath_1000G) 

export.bed(c(L1Ranges_shortName, L1NeighborRanges_1000G), L1NeighborbedPath_1000G) 

# Specify the general path to 1000 genome bam file
BamPath1000G_General <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/IndividualID/alignment/IndividualID.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam"
DirPath1000G_General <- "ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/IndividualID/alignment/"

# Get all files in the 1000 genome path
Files1000G <- list.files(Path1000G)
blnNotAnalyzed <- sapply(SampleColumns, function(x) length(grep(x, Files1000G)) == 0)
IndividualID <- SampleColumns[blnNotAnalyzed][3]
SampleColumns[blnNotAnalyzed]

########################################################################################
# Run MELT for low coverage bam files
IndividualID = SampleColumns[1]
for (IndividualID in SampleColumns[blnNotAnalyzed]){
  
  # Define paths
  DirPath1000G <- gsub("IndividualID", IndividualID, DirPath1000G_General)
  BamOutPath   <- paste(Path1000G, "L1Filtered_", IndividualID, ".bam", sep = "")
  
  # Create curl command to get file names on ftp directories
  DirListPath  <- paste(Path1000G, "FileList.txt", sep = "")
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
  BamPath1000G <- paste("ftp://", DirPath1000G, BamFileName, sep = "")

  # Construct samtools command to get filtered bam file
  SamToolsCmds <- c("module load samtools",
                    paste("samtools view -b -h", BamPath1000G, "-L", 
                          L1NeighborbedPath_1000G, "| samtools sort -o",
                          BamOutPath),
                    paste("samtools index", BamOutPath))
  
  # Construct MELT command
  MELTCmds <- paste("java -Xmx2g -jar /labs/dflev/hzudohna/MELTv2.1.5/MELT.jar Single -bamfile",
#                    BamOutPath,"-h /labs/dflev/hzudohna/RefSeqData/hg19masked.fa",
                    BamOutPath, 
                    "-c 7",
                    "-h /labs/dflev/hzudohna/RefSeqData/hs37d5.fa.gz",
                    "-t /labs/dflev/hzudohna/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip",
                    "-n /labs/dflev/hzudohna/MELTv2.1.5/add_bed_files/1KGP_Hg19/hg19.genes.bed ",
                    "-b MT/GL000207.1/GL000226.1/GL000229.1/GL000231.1/GL000210.1/GL000239.1/GL000235.1/GL000201.1/GL000247.1/GL000245.1/GL000197.1/GL000203.1/GL000246.1/GL000249.1/GL000196.1/GL000248.1/GL000244.1/GL000238.1/GL000202.1/GL000234.1/GL000232.1/GL000206.1/GL000240.1/GL000236.1/GL000241.1/GL000243.1/GL000242.1/GL000230.1/GL000237.1/GL000233.1/GL000204.1/GL000198.1/GL000208.1/GL000191.1/GL000227.1/GL000228.1/GL000214.1/GL000221.1/GL000209.1/GL000218.1/GL000220.1/GL000213.1/GL000211.1/GL000199.1/GL000217.1/GL000216.1/GL000215.1/GL000205.1/GL000219.1/GL000224.1/GL000223.1/GL000195.1/GL000212.1/GL000222.1/GL000200.1/GL000193.1/GL000194.1/GL000225.1/GL000192.1/NC_007605",
#                    "-bed", L1bedPath,
                    "-w /labs/dflev/hzudohna/1000Genomes/")
  
  # Command to remove bam files
  RemCmds <- c(paste("rm", BamOutPath), paste("rm ", BamOutPath, ".bai", sep = ""))
  
  ScriptName <- paste("ME_DEL_Script", IndividualID, sep = "_")
  CreateAndCallSlurmScript(file = ScriptName, 
                           SlurmCommandLines = c(SamToolsCmds, 
                                                 "echo 'bam file retrieved'", 
                                                 "module load bowtie",
                                                 MELTCmds, 
                                                 "echo 'MELT completed'",
                                                 RemCmds,
                                                 "echo 'bam files removed'"))
  
}

########################################################################################
# Run MELT for high coverage bam files

# Read file with high coverage paths and get IDs with 
HiCovPaths <- readLines("/labs/dflev/hzudohna/1000Genomes/1000G_highCoverPaths")
HighCoverageIDs <- sapply(HiCovPaths, function(x){
  PathSplit <- strsplit(x, "/")[[1]]
  strsplit(PathSplit[length(PathSplit)], "\\.")[[1]][1]
  
})
HighCoverageIDs
DirPath1000G_General_HiCov <- "ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/IndividualID/high_coverage_alignment/"
for (IndividualID in HighCoverageIDs){
  
  # Define paths
  DirPath1000G <- gsub("IndividualID", IndividualID, DirPath1000G_General_HiCov)
  BamOutPath   <- paste(Path1000G, "L1Filtered_HiCov_", IndividualID, ".bam", sep = "")
  
  # Create curl command to get file names on ftp directories
  DirListPath  <- paste(Path1000G, "FileList.txt", sep = "")
  CurlCmd <- paste("curl", paste("ftp://", DirPath1000G, sep = ""),
                   ">", DirListPath)
  system(CurlCmd)
  
  # Get name of bam file 
  DirListLines <- readLines(DirListPath)
  FileNames <- sapply(DirListLines, function(x){
    SplitLines <- strsplit(x, " ")[[1]]
    SplitLines[length(SplitLines)]
  })
  BamFileName <- FileNames[-grep("bam.b", FileNames)]
  BamPath1000G <- paste("ftp://", DirPath1000G, BamFileName, sep = "")
  
  # Construct samtools command to get filtered bam file
  SamToolsCmds <- c("module load samtools",
                    paste("samtools view -b -h", BamPath1000G, "-L", 
                          L1NeighborbedPath, "| samtools sort -o",
                          BamOutPath),
                    paste("samtools index", BamOutPath))
  
  # Construct MELT command
  MELTCmds <- paste("java -Xmx2g -jar /labs/dflev/hzudohna/MELTv2.1.5/MELT.jar Deletion-Genotype -bamfile",
                    #                    BamOutPath,"-h /labs/dflev/hzudohna/RefSeqData/hg19masked.fa",
                    BamOutPath, "-h /labs/dflev/hzudohna/RefSeqData/human_g1k_v37.fasta",
                    "-bed", L1bedPath,
                    "-w /labs/dflev/hzudohna/1000Genomes/")
  
  # Command to remove bam files
  RemCmds <- c(paste("rm", BamOutPath), paste("rm ", BamOutPath, ".bai", sep = ""))
  
  ScriptName <- paste("ME_DEL_Script", IndividualID, sep = "_")
  CreateAndCallSlurmScript(file = ScriptName, 
                           SlurmCommandLines = c(SamToolsCmds, 
                                                 "echo 'bam file retrieved'", 
                                                 MELTCmds, 
                                                 "echo 'MELT completed'",
                                                 RemCmds,
                                                 "echo 'bam files removed'"))
  
}

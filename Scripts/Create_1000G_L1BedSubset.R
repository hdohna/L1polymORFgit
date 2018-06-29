# The script below subsets creates a bed file to subset variants around L1s 
# from the 1000 Genome data
# Data were downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specify parameters
WindowWidth  <- 10^5
MinNrCarrier <- 5
NrInfoCol    <- 9
DataFolder  <- "/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/"
FilePattern <- "genotypes.vcf"
BedPath <- "/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/L1WindowSubset.bed"

# Read in 1000 genome L1 table
L1_1000G <- read.delim("/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/L1_1000G_withGenoNum",
                       header = T, sep = " ")

# Subset L1 file to retain only entries with enough carriers
blnL1 <- L1_1000G[,(NrInfoCol + 1):ncol(L1_1000G)] > 0
blnEnough <- rowSums(blnL1, na.rm = T) >= MinNrCarrier
sum(blnEnough)

# Create a bed file with windows around each L1
L1WindowsBed <- data.frame(chrom  = L1_1000G$CHROM[blnEnough],
                           chromStart = L1_1000G$POS[blnEnough] - WindowWidth,
                           chromEnd   = L1_1000G$POS[blnEnough] + WindowWidth
                           )

# Get index of overlapping rows
NR    <- nrow(L1WindowsBed)
idxOL <- which((L1WindowsBed$chromStart[-1] <= L1WindowsBed$chromEnd[-NR]) &
  (L1WindowsBed$chrom[-1] == L1WindowsBed$chrom[-NR]))

# Remove overlapping rows
L1WindowsBed$chromEnd[idxOL] <- L1WindowsBed$chromEnd[idxOL + 1]
L1WindowsBed <- L1WindowsBed[-(idxOL + 1), ]

# Write out bed file
write.table(L1WindowsBed, BedPath, quote = F, row.names = F)

# Get file names, loop over files and do the filtering
AllFiles <- list.files(DataFolder, pattern = FilePattern, full.names = T)
AllFiles <- AllFiles[-grep("vcf.", AllFiles)]
AllFiles

# Loop over file names, read file and append to existing
cat("Create L1 subset file per chromosome\n")
VcfFile <- AllFiles[1]

for (VcfFile in AllFiles[-1]){
  cat("Processing", VcfFile, "\n")
  Chrom <- strsplit(VcfFile, "\\.")[[1]][2]
  OutFile <- paste(DataFolder, Chrom, ".L1Windowsubset.vcf", sep = "")
  VcfCmd <- c("module load vcftools", 
              paste("vcftools --vcf", VcfFile, "--bed", BedPath,
              "--out", OutFile))
  ScriptName <- paste("L1Window", Chrom, sep = "_")
  CreateAndCallSlurmScript(file = ScriptName,
                           SlurmCommandLines = VcfCmd, 
                           scriptName = ScriptName) 
}

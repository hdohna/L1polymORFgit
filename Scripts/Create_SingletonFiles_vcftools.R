# The script below runs vcf commands to create files with singletons from the 1000 Genome data
# Data were downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specify parameters
NrInfoCols  <- 9
DataFolder  <- "/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/"
FilePattern <- "genotypes.vcf"

# Get file names, loop over files and do the filtering
AllFiles <- list.files(DataFolder, pattern = FilePattern, full.names = T)
AllFiles <- AllFiles[-grep("vcf.", AllFiles)]
AllFiles

# Loop over file names, read file and append to existing
cat("Create singleton file per chromosome\n")
VcfFile <- AllFiles[1]
for (VcfFile in AllFiles){
  cat("Processing", VcfFile, "\n")
  Chrom <- strsplit(VcfFile, "\\.")[[1]][2]
  VcfCmd <- c("module load vcftools", 
              paste("vcftools --vcf", VcfFile, "--singletons --out", Chrom))
  ScriptName <- paste("Singl", Chrom, sep = "_")
  CreateAndCallSlurmScript(file = ScriptName,
                           SlurmCommandLines = VcfCmd, 
                          scriptName = ScriptName) 
}


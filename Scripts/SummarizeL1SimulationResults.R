# The following script creates a summary of vcf files created by MELT on 
# simulated genomes with L1 insertions

load('/labs/dflev/hzudohna/RefSeqData/GRanges_L1_1000Genomes.RData')

# Specify simulation directory
SimDir <- "/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/"

# Get names of vcf files
VcfDirs <- list.dirs(SimDir, full.names = F)
VcfDirs <- VcfDirs[VcfDirs %in% SampleColumns]
VcfFiles <- paste(SimDir, VcfDirs, "/LINE1.final_comp.vcf", sep = "")

# Collect info on each file
FileInfo <- file.info(VcfFiles)
FileInfo$SampleID <- VcfDirs
write.csv(FileInfo, "/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/L1SimulationResults.csv")
# The following script adds read groups to all bam files in a folder

# Path to folder containing bam files
BamFolder <- "/home/hzudohna/L1polymORF/Data/PacbioFastqPerSuspectPeak/"

# Picard command
PicardCmd   <- "java -jar /home/txw/picard/picard-tools-1.52/picard.jar AddOrReplaceReadGroups"
OptionLines <- c("RGLB=lib1", "RGPL=illumina", "RGPU=unit1", "RGSM=20")

# Get all paths to fastq files in the folder
FileNames <- list.files(BamFolder, pattern = ".bam", full.names = T)
FileNames <- FileNames[-grep(".bam.", FileNames)]

# Create a command per file 
OptionLines <- paste(OptionLines, collapse = " ")
InFiles  <- paste("I=", FileNames, sep = "")
OutFiles <- paste("O=", FileNames, sep = "")
CmdLines <- paste(PicardCmd,  InFiles, OutFiles, OptionLines)

# Loop over command line and run them 
for (CmdL in CmdLines) system(CmdL)


SC50 <- SampleColumns[1:50]
SC50 [which(!SC50 %in% SampleIDs)]

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')


"HG00100" "HG00109" "HG00112" "HG00121" "HG00127" "HG00128" "HG00129"
"HG00133" "HG00138" "HG00141" "HG00142"

# Sam commands to turn sam file into bam file
SamBamCmds <- c("module load samtools",
                paste("samtools view /scratch/users/hzudohna/hg19_HG00109.bam -o", 
                      "/scratch/users/hzudohna/hg19_HG00109.sam"))


# Awk commands to replace MAPQ by 0 for unmapped reads
AwkCmds <- paste("awk '{if(int(($2 % 8) / 4) == 1){",
                "print $1, $2, $3, $4, 0, $6, $7, $8, $9, $10, $11;",
                "} else {",
                "print $0; } }'",
                "/scratch/users/hzudohna/hg19_HG00109.sam",
                "> /scratch/users/hzudohna/hg19_HG00109_sorted_MapQreplaced.sam")

# Sam commands to turn sam file into bam file
SamBamCmds2 <- c("module load samtools",
                paste("samtools view /scratch/users/hzudohna/hg19_HG00109.sam -b -o", 
                      "/scratch/users/hzudohna/hg19_HG00109_awk.bam"))

# Sam commands to sort and index bam file
SamIdxCmds <- c("module load samtools",
                paste("samtools sort /scratch/users/hzudohna/hg19_HG00109_awk.bam", 
                      "-o /scratch/users/hzudohna/hg19_HG00109_sorted_awk.bam"),
                paste("samtools index /scratch/users/hzudohna/hg19_HG00109_sorted_awk.bam"))

# Construct MELT command
MELTCmds <- c("module load bowtie2",
              paste("java -Xmx2g -jar /labs/dflev/hzudohna/MELTv2.1.5/MELT.jar Single -bamfile",
                    "/scratch/users/hzudohna/hg19_HG00109_sorted_MapQreplaced.sam", 
                    "-c 7",
                    "-h /labs/dflev/hzudohna/RefSeqData/hg19.fa", 
                    "-t /labs/dflev/hzudohna/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip",
                    "-n /labs/dflev/hzudohna/MELTv2.1.5/add_bed_files/1KGP_Hg19/hg19.genes.bed",
                    "-b MT/GL000207.1/GL000226.1/GL000229.1/GL000231.1/GL000210.1/GL000239.1/GL000235.1/GL000201.1/GL000247.1/GL000245.1/GL000197.1/GL000203.1/GL000246.1/GL000249.1/GL000196.1/GL000248.1/GL000244.1/GL000238.1/GL000202.1/GL000234.1/GL000232.1/GL000206.1/GL000240.1/GL000236.1/GL000241.1/GL000243.1/GL000242.1/GL000230.1/GL000237.1/GL000233.1/GL000204.1/GL000198.1/GL000208.1/GL000191.1/GL000227.1/GL000228.1/GL000214.1/GL000221.1/GL000209.1/GL000218.1/GL000220.1/GL000213.1/GL000211.1/GL000199.1/GL000217.1/GL000216.1/GL000215.1/GL000205.1/GL000219.1/GL000224.1/GL000223.1/GL000195.1/GL000212.1/GL000222.1/GL000200.1/GL000193.1/GL000194.1/GL000225.1/GL000192.1/NC_007605",
                    "-w /labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/HG00100"))


# Create a list of all commands
AllCmds <- c(SamBamCmds, SamIdxCmds, AwkCmds, MELTCmds)

# Create script name and run script
ScriptName <- "Sim_Awk_MELTScript"
ScriptName <- "Awk_Replace"
RunList <- CreateAndCallSlurmScript(file = ScriptName, 
                                    RunTime = '48:00:00',
                                    Mem     = '200G',
                                    SlurmCommandLines = AwkCmds)

RunList <- CreateAndCallSlurmScript(file = ScriptName, 
                                    RunTime = '48:00:00',
                                    Mem     = '200G',
                                    SlurmCommandLines = AllCmds)
cat(RunList$RunMessage, "\n")
RunIDs_SimMELT <- c(RunIDs_SimMELT, RunList$RunID)

# JobID: 10732996

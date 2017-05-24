# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

RL1 <- getReadLengthFromFastQ("/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/PacBioL1captPS/filtered_subreads.fastq")
RL2 <- getReadLengthFromFastQ("/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/FastqSubreads/BZ_NA12878L1capt5-9kb_sub_reads.fastq")
cat("******** saving results to /srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/PacBioReadLengths.RData *********\n")
save.image(file = "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/PacBioReadLengths.RData")


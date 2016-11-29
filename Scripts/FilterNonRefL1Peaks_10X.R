# The following script filters a 10X bam file to keep all reads that intersect
# with peaks away from reference L1
library(Rsamtools)

load("/srv/gsfs0/projects/levinson/hzudohna/10Xcapture/NA12878_L1Ranges10X.RData")

cat("Filtering NA12878_capt10X_aln2hg19.bam\n")
cat("Saving filtered file as NA12878_capt10X_aln2hg19_NonRefL1filtered.bam\n")
cat("in folder /srv/gsfs0/projects/levinson/hzudohna/10Xcapture/ \n")

param <- ScanBamParam(which = SuspectL1Ranges, what = scanBamWhat())
filterBam(file = "/srv/gsfs0/projects/levinson/hzudohna/10Xcapture/NA12878_capt10X_aln2hg19.bam", 
          destination = "/srv/gsfs0/projects/levinson/hzudohna/10Xcapture/NA12878_capt10X_aln2hg19_NonRefL1filtered.bam", 
          param = param)

java -XX:MaxHeapSize=1g -jar /srv/gsfs0/software/picard-tools/2.4.1/picard.jar AddOrReplaceReadGroups \
I=/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/CCS_BCDE.ccs.sorted.bam \
O=/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/CCS_BCDE.ccs.sorted.withRG.bam \
RGID=1 \
RGLB=lib1 \
RGPL=pacbio \
RGPU=unit1 \
RGSM=20
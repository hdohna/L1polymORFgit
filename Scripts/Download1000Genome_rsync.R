rsync -a -P rsync://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ ./
  
rsync -a -P rsync://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ ./srv/gsfs0/projects/levinson/hzudohna/1000Genomes/


ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz


rsync -P rsync://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ./srv/gsfs0/projects/levinson/hzudohna/1000Genomes/
rsync --copy-links --times --verbose rsync://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ./
#!/bin/bash
# Submit R batch job to queue
# Usage:
# Rbatch.sh [OPTIONS] infile [outfile]
#
# See R CMD BATCH --help for further details.
# By: Mario Pineda-Krch
#
# Note: the script has to be executable before it'll work on a unix box. Quite 
# likley you'll have to do 'chmod a+x Rbatch.sh' to set this after checking 
# out a working copy. 

ARGS="$*"
TPWD=$PWD
qsub -N $1 -V <<EOF
#!/bin/bash
#PBS -l nodes=1

# Join output log with error log
#PBS -j oe

cd $TPWD
pwd
echo "Running R with arguments:"
echo $ARGS
/share/diskarray1/Updated_RNA_software/R/R-3.2.1-download/bin/R CMD BATCH --no-save --no-restore $ARGS
EOF

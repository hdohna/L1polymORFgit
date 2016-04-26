# The script below runs the function 'AlignmentWorkflow'

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Run function
AlignmentWorkflow(FastqFile = "", 
                  ReferenceFile = "", 
                  SamFile = '/share/diskarray3/hzudohna/10XData/NA12878_10X_aln2L1.sam',
                  blnAlign = F,
                  blnSam2SortedBam = T,
                  blnDedup = T,
                  blnRemMapQ0 = T,
                  blnSort = T)
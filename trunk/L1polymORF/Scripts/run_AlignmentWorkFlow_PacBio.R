# The script below runs the function 'AlignmentWorkflow'

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Run function
AlignmentWorkflow(FastqFile = "", 
                  ReferenceFile = "", 
                  SamFile = '/share/diskarray3/hzudohna/NA12878PacBio_aln2Catalogue2016-05-07.sam',
                  blnAlign = F,
                  blnSam2SortedBam = T,
                  blnDedup = T,
                  blnRemMapQ0 = T,
                  blnSort = T)
# The script below runs the function 'AlignmentWorkflow'

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Run function
AlignmentWorkflow(FastqFile = "", 
                  ReferenceFile = "", 
                  SamFile = '/share/diskarray3/hzudohna/NA12878-L15P_S1_L001_001_Catalogue.sam',
                  BamFileSorted = NULL,
                  BamFileDedup = NULL,
                  blnAlign = F,
                  blnSam2SortedBam = T,
                  blnDedup = T,
                  blnSort = T)
# The script below runs the function 'AlignmentWorkflow'

# Collect command arguments (it is assumed to contain only the sam file name)
args <- commandArgs(trailingOnly = TRUE)

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Run function
cat("Running alignment workflow for", args[1])
AlignmentWorkflow(InFile = args[1], 
                  ReferenceFile = "", 
                  blnAlign = F,
                  blnSam2SortedBam = T,
                  blnDedup = T,
                  blnRemMapQ0 = T,
                  blnSort = T)
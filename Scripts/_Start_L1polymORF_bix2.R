##############################################
#
# General description:
#
#   The following script loads standard packages and functions to be used in
#   the L1polymORF project

# Input:
#
#   

# Output:
#   
#    : 

##############################################

########################################
#                                      #
#  Source packages and set parameters  #
#                                      #
########################################

# Source all functions from GeneralRFunctions folder
AllFunctions <- list.files(path = "/home/hzudohna/GeneralRFunctions/", 
                           pattern = ".[rR]", full.names = T)
sapply(AllFunctions, source)

# Source all functions from Functions folder
AllFunctions <- list.files(path = "/home/hzudohna/L1polymORF/Functions/", 
                           pattern = ".[rR]", full.names = T)
sapply(AllFunctions, source)

# Set paths to various NGS tools
system('export PATH=$PATH:/home/txw/samtools/samtools-1.2/')
system('export PATH=$PATH:/home/txw/bwa/bwa-0.7.12/')

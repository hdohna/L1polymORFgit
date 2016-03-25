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


# Load packages
library(seqinr)
library(ape)

# Source all functions from GeneralRFunctions folder
AllFunctions <- list.files(path = "D:/GeneralRFunctions/", 
                           pattern = ".[rR]", full.names = T)
sapply(AllFunctions, source)

# Source all functions from Functions folder
AllFunctions <- list.files(path = "D:/L1polymORF/Functions/", 
                           pattern = ".[rR]", full.names = T)
sapply(AllFunctions, source)
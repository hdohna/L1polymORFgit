##############################################
#
# General description:
#
#   The following script sources all R functions in the home directory
#   on SCG4 (Stanford genome cluster)


##############################################

########################################
#                                      #
#  Source all fu  #
#                                      #
########################################

# Source all functions from GeneralRFunctions folder
AllFunctions <- list.files(path = "/home/hzudohna/RFunctions/", 
                           pattern = ".[rR]", full.names = T)
sapply(AllFunctions, source)



##############################################
#
# General description:
#
#   The following function extracts the AIC from results produced by R function constrOptim
#   

# Input:
#
#     OptimResults: object produced by R function constrOptim


# Output:
#   
#    AIC value

# Comments:
#   

##############################################

GetAIC <- function(OptimResults){
  round(2 * (length(OptimResults$par) + OptimResults$value), 2)
}



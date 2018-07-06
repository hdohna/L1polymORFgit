##############################################
#
# General description:
#
#   The following function collapses a set of genomic ranges
#   

# Input:
#
#     Ranges: ranges to collapse
#     Width2Add: total width to add before collapsing


# Output:
#   
#    Collapsed granges

##############################################

CollapseGRanges <- function(Ranges, Width2Add){
  
  GR_Resized <- resize(Ranges, width = width(Ranges) + Width2Add, 
                       fix = "center")
  union(GR_Resized, GR_Resized)
}



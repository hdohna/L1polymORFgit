##############################################
#
# General description:
#
#   The following function filters a bam file (BamFileToBeFiltered) to retain
#   only reads that occur in file BamFilter

# Input:
#
#     BamFileToBeFiltered: path to bam file that is filtered
#     BamFilter: path to bam file that is used for filtering
#     OutFile: path to file that should be saved

# Output:
#   
#    Bam file saved under the name specified in OutFile

##############################################

######                                      
# Source packages and set parameters  
######

FilterBamByBam <- function(BamFileToBeFiltered, BamFilter){
  
  cat("Reading IDs of file", BamFilter,  "...\n")
  IDs <- scanBam(BamFilter, param=ScanBamParam(what="qname"))
  IDs <- unlist(IDs)
  
  cat("Filtering file", BamFileToBeFiltered, "by IDs of file", 
      BamFilter,  "\n")
  cat("Writing filtered file", OutFile, "\n")
  IDFilter <- FilterRules(getIDs <- function(DF){DF$qname %in% IDs})
  filterBam(BamFileToBeFiltered, OutFile, filter = IDFilter)
  
}


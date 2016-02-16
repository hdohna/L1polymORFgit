# General description:
#
#    This function takes a table and appends duplicated columns as rows

# Arguments:
#   
#    Table: table to be processed
#    MaxDupl: integer indicating maximum number of times columns might be
#       duplicated 

# Output:
#   
#    NewTable: table with duplicated columns appended as additional rows

# Comment:
#   
#   This function assumes that duplicated column names have an ".n" appended 
#   where "n" denotes the number of duplications (R does this when it reads in
#   tables)

AppendDuplColumns <- function(Table, MaxDupl = 5){
  
  # Create a StartTable with duplicated column names
  CNames        <- colnames(Table)
  DuplColNames  <- unlist(lapply(1:MaxDupl, function(i) paste(CNames, i, sep = ".")))
  blnNoDuplCols <- !CNames %in% DuplColNames
  NewTable      <- Table[,blnNoDuplCols]
  
  # Loop through duplicated column names and append them
  for (i in 1:MaxDupl) {
    DuplColNames  <- paste(CNames, i, sep = ".")
    blnDuplCols   <- CNames %in% DuplColNames
    TableAppend   <- Table[,blnDuplCols]
    CNamesDupl    <- CNames[blnDuplCols]
    CNamesNoDupl  <- substr(CNamesDupl, 1, nchar(CNamesDupl) - nchar(i) - 1)
    colnames(TableAppend) <- CNamesNoDupl
    NewTable      <- rbind(NewTable, TableAppend)
  }
  NewTable
}




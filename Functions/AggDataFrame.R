##############################################
#
# General description:
#
#   The following function aggregates a data frame
#   

# Input:
#
#     DF:        input data frame (required)
#     GroupCol:  name of grouping column (required)
#     MeanCols:  names of columns to be aggregated by mean (required)
#     SumCols:   names of columns to be summed
#     MedCols:   names of columns to be aggregated by median
#     LengthCols: names of columns to be aggregated by length
#     MaxCols: names of columns to be aggregated by maximum
#     MinCols: names of columns to be aggregated by minimum
#     VarCols: names of columns to be aggregated by length
#     Addcols:   name of additional columns to be added


# Output:
#   
#    Aggregated data frame

# Comments:
#   
#    

##############################################

AggDataFrame <- function(DF, GroupCol, MeanCols, 
                         SumCols = NULL,
                         MedCols = NULL, LengthCols = NULL,
                         MaxCols = NULL, MinCols = NULL,
                         VarCols = NULL, Addcols = NULL
){
  
  # Aggregate by mean (required)
  DF_Agg <- aggregate.data.frame(DF[, MeanCols], by = list(DF[, GroupCol]),
                                 FUN = function(x) mean(x, na.rm = T)) 
  colnames(DF_Agg) <- c(GroupCol, paste(MeanCols, "mean", sep = "_"))
  
  # Aggregate by sum (optional)
  if (!is.null(SumCols)) {
    DF_Sum <- aggregate.data.frame(DF[, SumCols], by = list(DF[, GroupCol]),
                                   FUN = function(x) sum(x, na.rm = T))
    colnames(DF_Sum) <- c(GroupCol, paste(SumCols, "sum", sep = "_"))
    DF_Agg <- merge(DF_Agg, DF_Sum)
  }

  # Aggregate by median (optional)
  if (!is.null(MedCols)) {
    DF_Med <- aggregate.data.frame(DF[, MedCols], by = list(DF[, GroupCol]),
                                   FUN = function(x) median(x, na.rm = T))
    colnames(DF_Med) <- c(GroupCol, paste(MedCols, "med", sep = "_"))
    DF_Agg <- merge(DF_Agg, DF_Med)
  }

  # Aggregate by length (optional)
  if (!is.null(LengthCols)) {
    DF_L <- aggregate.data.frame(DF[, LengthCols], by = list(DF[, GroupCol]),
                                   FUN = length)
    colnames(DF_L) <- c(GroupCol, paste(LengthCols, "N", sep = "_"))
    DF_Agg <- merge(DF_Agg, DF_L)
  }
  
  # Aggregate by maximum (optional)
  if (!is.null(MaxCols)) {
    DF_Max <- aggregate.data.frame(DF[, MaxCols], by = list(DF[, GroupCol]),
                                 FUN = function(x) max(x, na.rm = T))
    colnames(DF_Max) <- c(GroupCol, paste(MaxCols, "max", sep = "_"))
    DF_Agg <- merge(DF_Agg, DF_Max)
  }

  # Aggregate by minimum (optional)
  if (!is.null(MinCols)) {
    DF_Min <- aggregate.data.frame(DF[, MinCols], by = list(DF[, GroupCol]),
                                 FUN = function(x) min(x, na.rm = T))
    colnames(DF_Min) <- c(GroupCol, paste(MinCols, "min", sep = "_"))
    DF_Agg <- merge(DF_Agg, DF_Min)
  }

    # Aggregate by variance (optional)
  if (!is.null(VarCols)) {
    DF_Var <- aggregate.data.frame(DF[, VarCols], by = list(DF[, GroupCol]),
                                   FUN = function(x) var(x, na.rm = T))
    colnames(DF_Var) <- c(GroupCol, paste(VarCols, "var", sep = "_"))
    DF_Agg <- merge(DF_Agg, DF_Var)
  }
  
  # Add aditional columns
  if (!is.null(Addcols)) {
    NColBefore <- ncol(DF_Agg)
    idxMatch <- match(DF_Agg[ , GroupCol], DF[, GroupCol])
    DF_Agg   <- cbind(DF_Agg, DF[idxMatch, Addcols])
    colnames(DF_Agg)[(NColBefore + 1):(NColBefore + length(Addcols))] <- Addcols
  }
  DF_Agg
}



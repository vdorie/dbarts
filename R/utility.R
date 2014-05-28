## Turns data.frame w/factors into matrices of indicator variables. Dffers from
## model.matrix as it doesn't drop columns for co-linearity even with multiple
## factors. Also drops any non-numeric, non-factor columns.
makeModelMatrix <- function(x) {
  if (!is.data.frame(x)) warning("x in makeModelMatrix is not a dataframe");

  factorToMatrix <- function(factor) {
    sapply(levels(factor), { function(level) ifelse(factor == level, 1, 0) })
  }

  columnIsFactor <- sapply(x, is.factor)

  factorMatrices <- lapply((1:ncol(x))[columnIsFactor], function (colIndex) {
    result <- factorToMatrix(x[,colIndex])
    colnames(result) <- paste(colnames(x)[colIndex], ".", levels(x[,colIndex]), sep = "")
    result
  })
  factorNames <- unlist(lapply(factorMatrices, colnames))

  columnIsNumeric <- sapply(x, is.numeric)
  
  cbind(as.matrix(x[, columnIsNumeric]),
        matrix(unlist(factorMatrices), nrow(x), dimnames = list(NULL, factorNames)))
}

# Calculates the TI and WI index for networks. Input data frame should be without header, 
# the columns should be ordered: node name of species A, node name of species B, weight.
# Symmetrization_method is optional, can be either 'sum' (default), 'average', 'max', 'min' or 
# 'difference.' Multilink edges will be summed.

calculate_TI_WI <- function(input_data, steps, symmetrization_method="sum") {
  numstep <- steps
  symtype <- symmetrization_method
  
  if (!dim(input_data)[2] %in% c(2, 3)) stop("Wrong input format")
  
  if (is_character(input_data[1,1])) {
    names(input_data) = c("V1", "V2", "V3")
  }
  
  if (any(duplicated(input_data[,1:2]))) {
    input_data = aggregate(V3 ~ V1 + V2, data = input_data, FUN = sum)
    warning("Multilink edges were summed.")
  }
  
  nodeID <-
    levels(factor(c(
      unlist(input_data[, 1]), unlist(input_data[, 2])
    )))
  numnode <- length(nodeID)
  mx <- matrix(rep(0, numnode ^ 2), nrow = numnode, ncol = numnode)
  rownames(mx) <- nodeID
  colnames(mx) <- nodeID
  for (i in 1:length(input_data[, 1]))
  {
    mx[as.character(input_data[i, 1]), as.character(input_data[i, 2])] <- 1
    mx[as.character(input_data[i, 2]), as.character(input_data[i, 1])] <- 1
  }
  TI <- mx
  for (i in 1:numnode)
  {
    for (j in 1:numnode)
      TI[j, i] <- mx[j, i] / sum(mx[, i])
  }
  SI <- matrix(rep(0, numnode ^ 2), nrow = numnode, ncol = numnode)
  CI <- diag(numnode)
  for (i in 1:numstep)
  {
    CI <- CI %*% TI
    SI <- SI + CI
  }
  TI_index <- numeric(numnode)
  for (i in 1:numnode)
    TI_index[i] <- sum(SI[i,]) / numstep
  resu <- data.frame(nodeID, TI_index)
  
  
  if (dim(input_data)[2] == 2) {
    warning("No weights, computing only TI index.")
    return(resu)
  }
  
  if (!is.numeric(input_data[,3])) stop("Weights are not of numeric type")
  
  if (any(input_data[,3] < 0)) {
    input_data[,3] = abs(input_data[,3])
    warning("Negative values detected. Absolute values taken.")
  }
  
  mxw <- matrix(rep(0, numnode ^ 2), nrow = numnode, ncol = numnode)
  rownames(mxw) <- nodeID
  colnames(mxw) <- nodeID
  for (i in 1:length(input_data[, 1]))
  {
    mxw[as.character(input_data[i, 1]), as.character(input_data[i, 2])] <-
      input_data[i, 3]
  }
  for (i in 1:numnode)
    for (j in 1:numnode)
      if (i < j)
      {
        if ((mxw[i, j] == 0) ||
            (mxw[j, i] == 0)) {
          va <- mxw[i, j] + mxw[j, i]
          mxw[i, j] <- va
          mxw[j, i] <- va
        } else
        {
          if (symtype == "sum")
            va <- mxw[i, j] + mxw[j, i]
          if (symtype == "average")
            va <- (mxw[i, j] + mxw[j, i]) / 2
          if (symtype == "max")
            va <- max(mxw[i, j], mxw[j, i])
          if (symtype == "min")
            va <- min(mxw[i, j], mxw[j, i])
          if (symtype == "difference")
            va <- abs(mxw[i, j] - mxw[j, i])
          mxw[i, j] <- va
          mxw[j, i] <- va
        }
      }
  TI <- mxw
  for (i in 1:numnode)
  {
    for (j in 1:numnode)
      TI[j, i] <- mxw[j, i] / sum(mxw[, i])
  }
  SI <- matrix(rep(0, numnode ^ 2), nrow = numnode, ncol = numnode)
  CI <- diag(numnode)
  for (i in 1:numstep)
  {
    CI <- CI %*% TI
    SI <- SI + CI
  }
  WI_index <- numeric(numnode)
  for (i in 1:numnode)
    WI_index[i] <- sum(SI[i,]) / numstep
  resuw <- data.frame(nodeID, WI_index)
  
  result <- cbind(resu, WI_index)
  return(result)
}

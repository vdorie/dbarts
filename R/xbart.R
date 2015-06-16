xbart <- function(formula, data, subset, weights, offset, verbose = FALSE, n.samples = 200L,
                  K = 5L, n.reps = 200L, n.burn = c(100L, 20L), loss = c("mse", "mcr"), n.threads = guessNumCores(),
                  n.trees = 200L, k = 2, power = 2, base = 0.95, drop = TRUE,
                  resid.prior = chisq, control = dbartsControl(), sigma = NA_real_)
{
  matchedCall <- match.call()

  validateCall <- prepareCallWithArguments(matchedCall, quoteInNamespace(validateArgumentsInEnvironment), "control", "verbose", "n.samples", "sigma")
  validateCall <- addCallArgument(validateCall, 1L, sys.frame(sys.nframe()))
  eval(validateCall, parent.frame(1L), getNamespace("dbarts"))
  
  if (control@call != call("NA")[[1]]) control@call <- matchedCall
  control@verbose <- verbose
  
  dataCall <- prepareCallWithArguments(matchedCall, quoteInNamespace(dbartsData), "formula", "data", "subset", "weights", "offset")
  data <- eval(dataCall, parent.frame(1L))
  data@n.cuts <- rep_len(attr(control, "n.cuts"), ncol(data@x))
  data@sigma  <- sigma
  attr(control, "n.cuts") <- NULL
  
  if (is.na(data@sigma) && !control@binary)
    data@sigma <- summary(lm(data@y ~ data@x, weights = data@weights, offset = data@offset))$sigma
  
  
  uniqueResponses <- unique(data@y)
  if (length(uniqueResponses) == 2 && all(sort(uniqueResponses) == c(0, 1))) control@binary <- TRUE
  
  if (is.null(matchedCall$loss))
    loss <- loss[if (!control@binary) 1 else 2]
  
  tree.prior <- quote(cgm(power, base))
  tree.prior[[1]] <- quoteInNamespace(cgm)
  tree.prior[[2]] <- power[1]; tree.prior[[3]] <- base[1]
  tree.prior <- eval(tree.prior)

  node.prior <- quote(normal(k))
  node.prior[[1]] <- quoteInNamespace(normal)
  node.prior[[2]] <- k[1]
  node.prior <- eval(node.prior)
  
  resid.prior <-
    if (!is.null(matchedCall$resid.prior)) {
      eval(matchedCall$resid.prior, parent.frame(1L), getNamespace("dbarts"))
    } else {
      eval(formals(xbart)$resid.prior, getNamespace("dbarts"))()
    }
  model <- new("dbartsModel", tree.prior, node.prior, resid.prior)
  
  
  if (is.null(matchedCall$n.trees)) {
    n.trees <- control@n.trees
  } else {
    n.trees <- as.integer(n.trees)
    control@n.trees <- n.trees[1]
  }
  k       <- as.double(k)
  power   <- as.double(power)
  base    <- as.double(base)
  drop    <- as.logical(drop)
  
  result <- .Call(C_dbarts_xbart, control, model, data,
                  as.integer(K), as.integer(n.reps), as.integer(n.burn), loss, as.integer(n.threads),
                  n.trees, k, power, base, drop)
  
  if (is.null(result) || is.null(dim(result))) return(result)
  
  ## add dim names
  varNames <- c("n.trees", "k", "power", "base")
  if (identical(drop, TRUE))
    varNames <- varNames[sapply(varNames, function(varName) if (length(get(varName)) > 1) TRUE else FALSE)]
  
  dimNames <- vector("list", length(dim(result)))
  for (i in seq_along(varNames)) {
    x <- get(varNames[i])
    dimNames[[i]] <- as.character(if (is.double(x)) signif(x, 2) else x)
  }
  names(dimNames) <- c("rep", varNames)
  browser()
  dimnames(result) <- dimNames
  
  result
}

xbart <- function(formula, data, subset, weights, offset, verbose = FALSE, n.samples = 200L,
                  K = 5L, n.reps = 200L, n.burn = c(200L, 150L, 50L), loss = c("rmse", "mcr"),
                  n.threads = guessNumCores(),
                  n.trees = 200L, k = 2, power = 2, base = 0.95, drop = TRUE,
                  resid.prior = chisq, control = dbartsControl(), sigma = NA_real_)
{
  matchedCall <- match.call()

  validateCall <- prepareCallWithArguments(matchedCall, quoteInNamespace(validateArgumentsInEnvironment), "control", "verbose", "n.samples", "sigma")
  validateCall <- addCallArgument(validateCall, 1L, sys.frame(sys.nframe()))
  eval(validateCall, parent.frame(1L), getNamespace("dbarts"))
  
  if (control@call != call("NA")[[1L]]) control@call <- matchedCall
  control@verbose <- verbose
  
  dataCall <- prepareCallWithArguments(matchedCall, quoteInNamespace(dbartsData), "formula", "data", "subset", "weights", "offset")
  data <- eval(dataCall, parent.frame(1L))
  data@n.cuts <- rep_len(attr(control, "n.cuts"), ncol(data@x))
  data@sigma  <- sigma
  attr(control, "n.cuts") <- NULL
  
  kOrder <- order(k, decreasing = TRUE)
  kOrder.inv <- kOrder; kOrder.inv[kOrder] <- seq_along(kOrder)
  k <- k[kOrder]
  
  if (is.na(data@sigma) && !control@binary)
    data@sigma <- summary(lm(data@y ~ data@x, weights = data@weights, offset = data@offset))$sigma
  
  
  uniqueResponses <- unique(data@y)
  if (length(uniqueResponses) == 2L && all(sort(uniqueResponses) == c(0, 1))) control@binary <- TRUE
  
  if (is.null(matchedCall$loss)) {
    loss <- loss[if (!control@binary) 1L else 2L]
  } else if (is.function(loss)) {
    if (length(formals(loss)) != 2L) stop("supplied loss function must take exactly two arguments")
    loss <- list(loss, parent.frame(1L))
  } else if (is.list(loss)) {
    if (!is.function(loss[[1L]])) stop("first member of loss-list must be a function")
    if (length(formals(loss[[1L]])) != 2L) stop("supplied loss function must take exactly two arguments")
    if (!is.environment(loss[[2L]])) stop("second member of loss-list must be an environment")
  }
  
  tree.prior <- quote(cgm(power, base))
  tree.prior[[1L]] <- quoteInNamespace(cgm)
  tree.prior[[2L]] <- power[1L]; tree.prior[[3L]] <- base[1L]
  tree.prior <- eval(tree.prior)

  node.prior <- quote(normal(k))
  node.prior[[1L]] <- quoteInNamespace(normal)
  node.prior[[2L]] <- k[1L]
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
    n.trees <- coerceOrError(n.trees, "integer")
    control@n.trees <- n.trees[1L]
  }
  k       <- coerceOrError(k,     "numeric")
  power   <- coerceOrError(power, "numeric")
  base    <- coerceOrError(base,  "numeric")
  drop    <- coerceOrError(drop,  "logical")
  
  K         <- coerceOrError(K,         "integer")
  n.reps    <- coerceOrError(n.reps,    "integer")
  n.burn    <- coerceOrError(n.burn,    "integer")
  n.threads <- coerceOrError(n.threads, "integer")
  
  result <- .Call(C_dbarts_xbart, control, model, data,
                  K, n.reps, n.burn, loss,
                  n.threads, n.trees, k, power, base, drop)
  
  if (is.null(result) || is.null(dim(result))) return(result)
  
  ## add dim names
  varNames <- c("n.trees", "k", "power", "base")
  if (identical(drop, TRUE))
    varNames <- varNames[sapply(varNames, function(varName) if (length(get(varName)) > 1L) TRUE else FALSE)]
  
  if ("k" %in% varNames && any(kOrder != seq_along(k))) {
    indices <- rep(list(bquote()), length(dim(result)))
    indices[[1L + which(varNames == "k")]] <- kOrder.inv
    indexCall <- as.call(c(list(as.name("["), quote(result)), indices))
    result <- eval(indexCall)
    k <- k[kOrder.inv]
  }
  
  dimNames <- vector("list", length(dim(result)))
  for (i in seq_along(varNames)) {
    x <- get(varNames[i])
    dimNames[[i + 1L]] <- as.character(if (is.double(x)) signif(x, 2L) else x)
  }
  names(dimNames) <- c("rep", varNames)
  dimnames(result) <- dimNames
  
  result
}

xbart <- function(formula, data, subset, weights, offset, verbose = FALSE, n.samples = 200L,
                  method = c("k-fold", "random subsample"), n.test = c(5, 0.2),
                  n.reps = 40L, n.burn = c(200L, 150L, 50L), loss = c("rmse", "log", "mcr"),
                  n.threads = guessNumCores(),
                  n.trees = 75L, k = 2, power = 2, base = 0.95, drop = TRUE,
                  resid.prior = chisq, control = dbartsControl(), sigma = NA_real_)
{
  matchedCall <- match.call()
  
  currEnv <- sys.frame(sys.nframe())
  evalEnv <- parent.frame(1L)
  
  validateCall <- redirectCall(matchedCall, quoteInNamespace(validateArgumentsInEnvironment), control, verbose, n.samples, sigma)
  validateCall <- addCallArgument(validateCall, 1L, currEnv)
  validateCall <- addCallArgument(validateCall, 2L, xbart)
  eval(validateCall, evalEnv, getNamespace("dbarts"))
  
  if (control@call != call("NA")[[1L]]) control@call <- matchedCall
  control@verbose <- verbose
  control@n.chains <- 1L
  control@keepTrees <- FALSE
  
  dataCall <- redirectCall(matchedCall, quoteInNamespace(dbartsData), formula, data, subset, weights, offset)
  data <- eval(dataCall, evalEnv)
  data@n.cuts <- rep_len(attr(control, "n.cuts"), ncol(data@x))
  data@sigma  <- sigma
  attr(control, "n.cuts") <- NULL
  
  uniqueResponses <- unique(data@y)
  if (length(uniqueResponses) == 2L && all(sort(uniqueResponses) == c(0, 1))) control@binary <- TRUE
  
  if (is.na(data@sigma) && !control@binary)
    data@sigma <- summary(lm(data@y ~ data@x, weights = data@weights, offset = data@offset))$sigma
  
  if (control@binary && is.null(matchedCall[["resid.prior"]]))
    matchedCall[["resid.prior"]] <- quote(fixed(1))
  
  if (!is.character(method) || method[1L] %not_in% eval(formals(xbart)$method))
    stop("method must be in '", paste0(eval(formals(xbart)$method), collapse = "', '"), "'")
  method <- method[1L]
  if (!is.null(matchedCall$method) && is.null(matchedCall$n.test))
    n.test <- eval(formals(xbart)$n.test)[match(method, eval(formals(xbart)$method))]
  n.test <- n.test[1L]
  
  if (is.null(matchedCall$loss)) {
    loss <- loss[if (!control@binary) 1L else 2L]
  } else if (is.function(loss)) {
    if (length(formals(loss)) != 2L) stop("supplied loss function must take exactly two arguments")
    loss <- list(loss, evalEnv)
  } else if (is.list(loss)) {
    if (!is.function(loss[[1L]])) stop("first member of loss-list must be a function")
    if (length(formals(loss[[1L]])) != 2L) stop("supplied loss function must take exactly two arguments")
    if (!is.environment(loss[[2L]])) stop("second member of loss-list must be an environment")
  }
    
  if (is.null(matchedCall$n.trees) && "n.trees" %not_in% names(matchedCall)) {
    n.trees <- control@n.trees
  } else {
    n.trees <- coerceOrError(n.trees, "integer")
    control@n.trees <- n.trees[1L]
  }
  
  k      <- coerceOrError(k, "numeric")
  kOrder <- order(k, decreasing = TRUE)
  kOrder.inv <- kOrder; kOrder.inv[kOrder] <- seq_along(kOrder)
  k <- k[kOrder]
  
  power   <- coerceOrError(power, "numeric")
  base    <- coerceOrError(base,  "numeric")
  drop    <- coerceOrError(drop,  "logical")
  
  tree.prior <- quote(cgm(power, base))
  tree.prior[[1L]] <- quoteInNamespace(cgm)
  tree.prior[[2L]] <- power[1L]; tree.prior[[3L]] <- base[1L]
  tree.prior <- eval(tree.prior)

  node.prior <- quote(normal(k))
  node.prior[[1L]] <- quoteInNamespace(normal)
  node.prior[[2L]] <- k[1L]
  node.prior <- eval(node.prior)
  
  resid.prior <-
    if (!is.null(matchedCall$resid.prior) || "resid.prior" %in% names(matchedCall)) {
      eval(matchedCall$resid.prior, evalEnv, getNamespace("dbarts"))
    } else {
      eval(formals(xbart)$resid.prior, getNamespace("dbarts"))()
    }
  model <- new("dbartsModel", tree.prior, node.prior, resid.prior,
               node.scale = if (control@binary) 3.0 else 0.5)
  
  if (method == "k-fold") {
    n.test <- coerceOrError(n.test, "integer")
    if (n.test < 2L || n.test > length(data@y))
      stop("for k-fold crossvalidation, n.test must be an integer in [2, N]")
  } else {
    n.test <- coerceOrError(n.test, "numeric")
    if (n.test > 1) n.test <- n.test / length(data@y)
    if (n.test <= 0 || n.test >= 1) stop("for random subsample crossvalidation, n.test must be in (0, 1)")
  }
  
  n.reps    <- coerceOrError(n.reps,    "integer")
  n.burn    <- coerceOrError(n.burn,    "integer")
  n.threads <- coerceOrError(n.threads, "integer")
  
  result <- .Call(C_dbarts_xbart, control, model, data, method,
                  n.test, n.reps, n.burn, loss,
                  n.threads, n.trees, k, power, base, drop)
  
  if (is.null(result) || is.null(dim(result))) return(result)
  
  ## add dim names
  varNames <- c("n.trees", "k", "power", "base")
  if (identical(drop, TRUE))
    varNames <- varNames[sapply(varNames, function(varName) if (length(get(varName)) > 1L) TRUE else FALSE)]
  
  if ("k" %in% varNames && any(kOrder != seq_along(k))) {
    indices <- rep(list(bquote()), length(dim(result)))
    indices[[1L + which(varNames == "k")]] <- kOrder.inv
    indexCall <- as.call(c(list(as.name("["), quote(result)), indices, drop = FALSE))
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

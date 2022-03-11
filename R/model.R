setMethod("initialize", "dbartsModel",
          function(.Object, tree.prior, node.prior, node.hyperprior, resid.prior,
                   proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4,
                                      birth = 0.5),
                   node.scale = 0.5)
{
  if (!missing(tree.prior)) .Object@tree.prior  <- tree.prior
  if (!missing(node.prior)) .Object@node.prior  <- node.prior
  if (!missing(node.hyperprior)) .Object@node.hyperprior  <- node.hyperprior
  if (!missing(resid.prior)) .Object@resid.prior <- resid.prior
  
  if (is.null(proposal.probs))
    proposal.probs <- c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)
  
  probs <- proposal.probs[c("birth_death", "swap", "change")]
  if (sum(is.na(probs)) == 1L) {
    probs[is.na(probs)] <- 1 - sum(probs[!is.na(probs)])
    names(probs) <- c("birth_death", "swap", "change")
  } else if (all(is.na(probs))) {
    probs <- c(birth_death = 0.5, swap = 0.1, change = 0.4)
  }
  
  .Object@p.birth_death <- probs[["birth_death"]]
  .Object@p.swap <- probs[["swap"]]
  .Object@p.change <- probs[["change"]]
  
  probs <- proposal.probs["birth"]
  if (is.na(probs)) probs <- c(birth = 0.5)
  
  
  .Object@p.birth <- probs[["birth"]]
  
  .Object@node.scale <- node.scale

  validObject(.Object)
  .Object
})

parsePriors <- function(control, data, tree.prior, node.prior, resid.prior, parentEnv)
{
  matchedCall <- match.call()

  evalEnv <- new.env(parent = parentEnv)
  evalEnv$control <- control
  evalEnv$data <- data
  evalEnv$cgm <- cgm
  evalEnv$normal <- normal
  evalEnv$chisq <- chisq
  evalEnv$fixed <- fixed
  
  # sub in a different default for node prior if data are binary
  if (control@binary)
    formals(evalEnv$normal)[["k"]] <- quote(chi(1.25, Inf))

  if (is.symbol(matchedCall$tree.prior))
    matchedCall$tree.prior  <- call(as.character(matchedCall$tree.prior))
  if (is.symbol(matchedCall$resid.prior))
    matchedCall$resid.prior <- call(as.character(matchedCall$resid.prior))
  if (is.symbol(matchedCall$node.prior))
    matchedCall$node.prior  <- call(as.character(matchedCall$node.prior))
  
  tree.prior  <- eval(matchedCall$tree.prior, evalEnv)
  resid.prior <- eval(matchedCall$resid.prior, evalEnv)
  
  node.prior <- node.hyperprior <- NULL
  massign[node.prior, node.hyperprior] <- eval(matchedCall$node.prior, evalEnv)
  
  namedList(tree.prior, resid.prior, node.prior, node.hyperprior)
}


cgm <- function(power = 2, base = 0.95, split.probs = 1 / num.vars)
{
  matchedCall <- match.call()
  if (is.null(matchedCall$split.probs))
    matchedCall$split.probs <- formals(quoteInNamespace(cgm))$split.probs
  split.probs.expr <- subTermInLanguage(matchedCall$split.probs, quote(num.vars), quote(ncol(data@x)))
  split.probs.expr <- subTermInLanguage(split.probs.expr, quote(numvars), quote(ncol(data@x)))
  split.probs <- eval(split.probs.expr, parent.frame())
  data <- parent.frame()$data
  
  if (length(split.probs) == 1L) {
    # if length 1, we can ignore it
    split.probs <- numeric()
  } else if (!is.null(names(split.probs))) {
    default <- NA_real_
    split.names <- names(split.probs)
    defaultMatch <- split.names %in% ".default"
    if (sum(defaultMatch) > 1L)
      stop("cannot assign split probabilities: default specified multiple times")
    if (sum(defaultMatch) == 1L) {
      default <- split.probs[[which(defaultMatch)]]
      split.probs <- split.probs[!defaultMatch]
      split.names <- names(split.probs)
    }
    
    result <- rep(default, ncol(data@x))
    names(result) <- colnames(data@x)

    if (is.null(names(result)) && length(split.names) > 0L)
      stop("cannot assign split probabilities: model matrix has no column names")
    
    namesMatch <- match(split.names, names(result))
    result[namesMatch[!is.na(namesMatch)]] <- split.probs[!is.na(namesMatch)]

    split.probs <- split.probs[is.na(namesMatch)]
    split.names <- names(split.probs)
    
    for (i in seq_along(split.probs)) {
      if (split.names[i] %not_in% attr(data@x, "term.labels"))
        stop("cannot assign split probabilities: unrecognized variable name '", split.names[i], "'")
      factorMatch <- which(startsWith(names(result), paste0(split.names[i], ".")))
      result[factorMatch] <- split.probs[i]
    }

    split.probs <- result

  } else {
    if (length(split.probs) != ncol(data@x))
      stop("cannot assign split probabilities: length of input (", length(split.probs),
           ") does not equal number of columns in model matrix (", ncol(data@x), ")")
  }
  
  if (length(split.probs) > 0L)
    split.probs <- split.probs / sum(split.probs)

  new("dbartsCGMPrior",
      power = power,
      base = base,
      splitProbabilities = split.probs)
}

normal <- function(k = 2.0)
{
  matchedCall <- match.call()
  evalEnv <- new.env(parent = parent.frame())
  evalEnv$chi <- function(degreesOfFreedom = 1.25, scale = Inf)
    new("dbartsChiHyperprior", degreesOfFreedom = degreesOfFreedom, scale = scale)

  if (!is.null(matchedCall[["k"]])) {
    kExpr <- matchedCall[["k"]]
    for (i in seq_len(2L)) {
      if (is.numeric(kExpr) || inherits(kExpr, "dbartsNodeHyperprior")) break
      
      if (is.character(kExpr)) {
        if (startsWith(kExpr, "chi")) {
          kExpr <- parse(text = kExpr)[[1L]]
          if (!is.call(kExpr))
            kExpr <- call(as.character(kExpr))
        }
        else kExpr <- coerceOrError(kExpr, "numeric")
      }
      if (is.symbol(kExpr) && !is.call(kExpr) && startsWith(as.character(kExpr), "chi"))
        kExpr <- call(as.character(kExpr))
      
      # the below evaluation might only lead to a lookup, in which case we have to do an
      # additional level of casting/eval
      kExpr <- eval(kExpr, evalEnv)
    }
    kHyperprior <- kExpr
  } else {
    kHyperprior <- eval(formals()[["k"]], evalEnv)
  }
  
  if (is.numeric(kHyperprior))
    kHyperprior <- new("dbartsFixedHyperprior", k = kHyperprior)
  
  namedList(node.prior = new("dbartsNormalPrior"), node.hyperprior = kHyperprior)
}

chisq <- function(df = 3, quant = 0.9)
{
  new("dbartsChiSqPrior", df = df, quantile = quant)
}
    
fixed <- function(value = 1.0)
{
  new("dbartsFixedPrior", value = value)
}


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


cgm <- function(power = 2, base = 0.95)
{
  new("dbartsCGMPrior", power = power, base = base)
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


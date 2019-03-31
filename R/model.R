setMethod("initialize", "dbartsModel",
          function(.Object, tree.prior, node.prior, resid.prior,
                   p.birth_death = 0.5, p.swap = 0.1, p.change = 0.4,
                   p.birth = 0.5, node.scale = 0.5)
{
  if (!missing(tree.prior)) .Object@tree.prior  <- tree.prior
  if (!missing(node.prior)) .Object@node.prior  <- node.prior
  if (!missing(resid.prior)) .Object@resid.prior <- resid.prior
  
  .Object@p.birth_death <- p.birth_death
  .Object@p.swap <- p.swap
  .Object@p.change <- p.change
  .Object@p.birth <- p.birth
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

  if (is.symbol(matchedCall$tree.prior)) matchedCall$tree.prior <- call(as.character(matchedCall$tree.prior))
  if (is.symbol(matchedCall$resid.prior)) matchedCall$resid.prior <- call(as.character(matchedCall$resid.prior))
  if (is.symbol(matchedCall$node.prior)) matchedCall$node.prior <- call(as.character(matchedCall$node.prior))
  
  tree.prior <- eval(matchedCall$tree.prior, evalEnv)
  resid.prior <- eval(matchedCall$resid.prior, evalEnv)
  node.prior <- eval(matchedCall$node.prior, evalEnv)
  
  namedList(tree.prior, resid.prior, node.prior)
}


cgm <- function(power = 2, base = 0.95)
{
  new("dbartsCGMPrior", power = power, base = base)
}

normal <- function(k = 2.0)
{
  matchedCall <- match.call()
  
  if (!is.null(matchedCall[["k"]])) {
    evalEnv <- new.env(parent = parent.frame())
    evalEnv$chi <- function(degreesOfFreedom = 2.5, scale = 1)
      new("dbartsChiHyperprior", degreesOfFreedom = degreesOfFreedom, scale = scale)
    
    kExpr <- matchedCall[["k"]]
    if (is.symbol(kExpr) && startsWith(as.character(kExpr), "chi")) {
      kExpr <- call(as.character(kExpr))
    } else if (is.character(kExpr) && startsWith(kExpr, "chi")) {
      kExpr <- parse(text = kExpr)[[1L]]
      if (!is.call(kExpr))
        kExpr <- call(as.character(kExpr))
    }
    
    k <- eval(kExpr, evalEnv)
  }
  
  new("dbartsNormalPrior", k = k)
}

chisq <- function(df = 3, quant = 0.9)
{
  new("dbartsChiSqPrior", df = df, quantile = quant)
}
    
fixed <- function(value = 1.0)
{
  new("dbartsFixedPrior", value = value)
}


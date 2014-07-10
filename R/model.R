setMethod("initialize", "dbartsModel",
          function(.Object, tree.prior, node.prior, resid.prior,
                   p.birth_death = 0.5, p.swap = 0.1, p.change = 0.4,
                   p.birth = 0.5)
{
  if (!missing(tree.prior)) .Object@tree.prior  <- tree.prior
  if (!missing(node.prior)) .Object@node.prior  <- node.prior
  if (!missing(resid.prior)) .Object@resid.prior <- resid.prior
  
  .Object@p.birth_death <- p.birth_death
  .Object@p.swap <- p.swap
  .Object@p.change <- p.change
  .Object@p.birth <- p.birth

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

  if (is.symbol(matchedCall$tree.prior)) matchedCall$tree.prior <- call(as.character(matchedCall$tree.prior))
  if (is.symbol(matchedCall$node.prior)) matchedCall$node.prior <- call(as.character(matchedCall$node.prior))
  if (is.symbol(matchedCall$resid.prior)) matchedCall$resid.prior <- call(as.character(matchedCall$resid.prior))
  
  tree.prior  <- eval(matchedCall$tree.prior, evalEnv)
  node.prior  <- eval(matchedCall$node.prior, evalEnv)
  resid.prior <- eval(matchedCall$resid.prior, evalEnv)
  
  namedList(tree.prior, node.prior, resid.prior)
}


cgm <- function(power = 2, base = 0.95)
{
  new("dbartsCGMPrior", power = power, base = base)
}

normal <- function(k = 2.0)
{
  new("dbartsNormalPrior", k = k)
}

chisq <- function(df = 3, quant = 0.9)
{
  new("dbartsChiSqPrior", df = df, quantile = quant)
}
    



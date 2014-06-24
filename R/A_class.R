## stupid file name to load first

setClass("cbartTreePrior")
setClass("cbartCGMPrior", contains = "cbartTreePrior",
         slots = list(power = "numeric", base = "numeric"),
         validity = function(object) {
           if (object@power <= 0.0) return("'power' must be positive.")
           if (object@base <= 0.0 || object@base >= 1.0) return("'base' must be in (0, 1).")
         })

setClass("cbartNodePrior")
setClass("cbartNormalPrior", contains = "cbartNodePrior",
         slots = list(k = "numeric"),
         validity = function(object) {
           if (object@k <= 0.0) return("'k' must be positive.")
         })

setClass("cbartResidPrior")
setClass("cbartChiSqPrior", contains = "cbartResidPrior",
         slots = list(df = "numeric", quantile = "numeric"),
         validity = function(object) {
           if (object@df <= 0.0) return("'df' must be positive.")
           if (object@quantile <= 0.0) return("'quantile' must be positive.")
         })

cbartControl <-
  setClass("cbartControl",
           slots =
           list(binary           = "logical",
                verbose          = "logical",
                keepTrainingFits = "logical",
                useQuantiles     = "logical",
                n.samples        = "integer",
                n.burn           = "integer",
                n.trees          = "integer",
                n.threads        = "integer",
                n.thin           = "integer",
                printEvery       = "integer",
                printCutoffs     = "integer",
                call             = "language"),
           prototype =
           list(binary           = FALSE,
                verbose          = TRUE,
                keepTrainingFits = TRUE,
                useQuantiles     = FALSE,
                n.samples        = NA_integer_,
                n.burn           = 100L,
                n.trees          = 200L,
                n.threads        = 1L,
                n.thin           = 1L,
                printEvery       = 100L,
                printCutoffs     = 0L,
                call             = quote(quote(NA()))),
           validity = function(object) {
             if (is.na(object@keepTrainingFits)) return("'keepTrainingFits' must be TRUE/FALSE.")
             if (is.na(object@useQuantiles))     return("'useQuantiles' must be TRUE/FALSE.")
             
             if (object@n.burn    <  0L) return("'n.burn' must be a non-negative integer.")
             if (object@n.trees   <= 0L) return("'n.trees' must be a positive integer.")
             if (object@n.threads <= 0L) return("'n.threads' must be a positive integer.")
             if (object@n.thin    <  0L) return("'n.thin' must be a non-negative integer.")
             
             if (object@printEvery   < 0L) return("'printEvery' must be a non-negative integer.")
             if (object@printCutoffs < 0L) return("'printCutoffs' must be a non-negative integer.")

             ## handle this in particular b/c it is set through cbart, not
             ## standard initializer
             if (!is.na(object@n.samples) && object@n.samples <= 0L) return("'n.samples' must be a positive integer.")

             TRUE
           })

setClass("cbartModel",
         slots =
         list(p.birth_death = "numeric",
              p.swap        = "numeric",
              p.change      = "numeric",
              
              p.birth = "numeric",

              tree.prior  = "cbartTreePrior",
              node.prior  = "cbartNodePrior",
              resid.prior = "cbartResidPrior"),
         prototype =
         list(p.birth_death = 1.0,
              p.swap        = 0.0,
              p.change      = 0.0,
              p.birth       = 0.5,
              tree.prior    = new("cbartCGMPrior"),
              node.prior    = new("cbartNormalPrior"),
              resid.prior   = new("cbartChiSqPrior")),
         validity = function(object) {
           proposalProbs <- c(object@p.birth_death, object@p.swap, object@p.change)
           if (any(proposalProbs < 0.0) || any(proposalProbs > 1.0)) return("Rule proposal probabilities must be in [0, 1].")
           if (abs(sum(proposalProbs) - 1.0) >= 1e-10) return("Rule proposal probabilities must sum to 1")
           
           if (object@p.birth <= 0.0 || object@p.birth >= 1.0) return("Birth probability for birth/death step must be in (0, 1).")
         })

setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("numericOrNULL", c("numeric", "NULL"))

setClass("cbartData",
         slots =
         list(y        = "numeric",
              x        = "matrix",
              varTypes = "integer",
              x.test   = "matrixOrNULL",
              weights  = "numericOrNULL",
              offset   = "numericOrNULL",
              n.cuts   = "integer",
              sigma    = "numeric"),
         prototype =
         list(y        = numeric(0),
              x        = matrix(0, 0, 0),
              varTypes = integer(0),
              x.test   = NULL,
              weights  = NULL,
              offset   = NULL,
              n.cuts   = integer(0),
              sigma    = NA_real_),
         validity = function(object) {
           numObservations <- length(object@y)
           if (nrow(object@x) != numObservations) return("Number of rows of 'x' must equal length of 'y'.")
           
           if (length(object@varTypes) > 0 &&
               any(object@varTypes != ORDINAL_VARIABLE && object@varTypes != CATEGORICAL_VARIABLE))
             return("Variable types must all be ordinal or categorical.")
           
           if (!is.null(object@x.test) && ncol(object@x.test) != ncol(object@x)) return("'x.test' must be null or have number of columns equal to 'x'.")
           if (!is.null(object@weights)) {
             if (length(object@weights) != numObservations) return("'weights' must be null or have length equal to that of 'y'.")
             if (anyNA(object@weights)) return("'weights' cannot be NA.")
             if (any(object@weights <= 0.0)) return("'weights' must all be positive.")
           }
           if (!is.null(object@offset) && length(object@offset) != numObservations) return("'offset' must be null or have length equal to that of 'y'.")
           if (!is.na(object@n.cuts) && length(object@n.cuts) != ncol(object@x)) return("Length of 'n.cuts' must equal number of columns in 'x'.")
           
           if (!is.na(object@sigma) && object@sigma <= 0.0) return("'sigma' must be positive.")
         })

## this shouldn't ever get created, used, modified, whathaveyou
setClass("cbartState",
         slots =
         list(fit.tree  = "numeric",
              fit.total = "numeric",
              fit.test  = "numeric",
              sigma     = "numeric",
              trees     = "character"))
              

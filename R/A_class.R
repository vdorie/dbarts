## stupid file name to load first

setClass("dbartsTreePrior")
setClass("dbartsCGMPrior", contains = "dbartsTreePrior",
         slots = list(power = "numeric", base = "numeric"),
         validity = function(object) {
           if (object@power <= 0.0) return("'power' must be positive")
           if (object@base  <= 0.0 || object@base >= 1.0) return("'base' must be in (0, 1)")
         })

setClass("dbartsNodePrior")
setClass("dbartsNormalPrior", contains = "dbartsNodePrior",
         slots = list(k = "numeric"),
         validity = function(object) {
           if (object@k <= 0.0) return("'k' must be positive")
         })

setClass("dbartsResidPrior")
setClass("dbartsChiSqPrior", contains = "dbartsResidPrior",
         slots = list(df = "numeric", quantile = "numeric"),
         validity = function(object) {
           if (object@df <= 0.0) return("'df' must be positive")
           if (object@quantile <= 0.0) return("'quantile' must be positive")
         })

setClass("dbartsControl",
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
              updateState      = "logical",
              call             = "language"),
         prototype =
         list(binary           = FALSE,
              verbose          = FALSE,
              keepTrainingFits = TRUE,
              useQuantiles     = FALSE,
              n.samples        = NA_integer_,
              n.burn           = 100L,
              n.trees          = 200L,
              n.threads        = 1L,
              n.thin           = 1L,
              printEvery       = 100L,
              printCutoffs     = 0L,
              updateState      = TRUE,
              call = quote(call("NA"))),
##              call             = quote(quote(NA()))),
         validity = function(object) {
           if (is.na(object@verbose))          return("'verbose' must be TRUE/FALSE")
           if (is.na(object@keepTrainingFits)) return("'keepTrainingFits' must be TRUE/FALSE")
           if (is.na(object@useQuantiles))     return("'useQuantiles' must be TRUE/FALSE")
           
           if (object@n.burn    <  0L) return("'n.burn' must be a non-negative integer")
           if (object@n.trees   <= 0L) return("'n.trees' must be a positive integer")
           if (object@n.threads <= 0L) return("'n.threads' must be a positive integer")
           if (object@n.thin    <  0L) return("'n.thin' must be a non-negative integer")
           
           if (object@printEvery   < 0L) return("'printEvery' must be a non-negative integer")
           if (object@printCutoffs < 0L) return("'printCutoffs' must be a non-negative integer")
           
           if (is.na(object@updateState)) return("'updateState' must be TRUE/FALSE")
           
           ## handle this in particular b/c it is set through dbarts, not
           ## standard initializer
           if (!is.na(object@n.samples) && object@n.samples < 0L) return("'n.samples' must be a non-negative integer")
           
           TRUE
         })

setClass("dbartsModel",
         slots =
         list(p.birth_death = "numeric",
              p.swap        = "numeric",
              p.change      = "numeric",
              
              p.birth = "numeric",

              tree.prior  = "dbartsTreePrior",
              node.prior  = "dbartsNodePrior",
              resid.prior = "dbartsResidPrior"),
         prototype =
         list(p.birth_death = 1.0,
              p.swap        = 0.0,
              p.change      = 0.0,
              p.birth       = 0.5,
              tree.prior    = new("dbartsCGMPrior"),
              node.prior    = new("dbartsNormalPrior"),
              resid.prior   = new("dbartsChiSqPrior")),
         validity = function(object) {
           proposalProbs <- c(object@p.birth_death, object@p.swap, object@p.change)
           if (any(proposalProbs < 0.0) || any(proposalProbs > 1.0)) return("rule proposal probabilities must be in [0, 1]")
           if (abs(sum(proposalProbs) - 1.0) >= 1e-10) return("rule proposal probabilities must sum to 1")
           
           if (object@p.birth <= 0.0 || object@p.birth >= 1.0) return("birth probability for birth/death step must be in (0, 1)")
         })

setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("numericOrNULL", c("numeric", "NULL"))

setClass("dbartsData",
         slots =
         list(y           = "numeric",
              x           = "matrix",
              varTypes    = "integer",
              x.test      = "matrixOrNULL",
              weights     = "numericOrNULL",
              offset      = "numericOrNULL",
              offset.test = "numericOrNULL",
              n.cuts      = "integer",
              sigma       = "numeric",

              testUsesRegularOffset   = "logical"),
         prototype =
         list(y           = numeric(0),
              x           = matrix(0, 0, 0),
              varTypes    = integer(0),
              x.test      = NULL,
              weights     = NULL,
              offset      = NULL,
              offset.test = NULL,
              n.cuts      = integer(0),
              sigma       = NA_real_,

              testUsesRegularOffset   = NA
              ),
         validity = function(object) {
           numObservations <- length(object@y)
           if (nrow(object@x) != numObservations) return("number of rows of 'x' must equal length of 'y'")
           
           if (length(object@varTypes) > 0 &&
               any(object@varTypes != ORDINAL_VARIABLE && object@varTypes != CATEGORICAL_VARIABLE))
             return("variable types must all be ordinal or categorical")
           
           if (!is.null(object@x.test) && ncol(object@x.test) != ncol(object@x)) return("'x.test' must be null or have number of columns equal to 'x'")
           if (!is.null(object@weights)) {
             if (length(object@weights) != numObservations) return("'weights' must be null or have length equal to that of 'y'")
             if (anyNA(object@weights)) return("'weights' cannot be NA")
             if (any(object@weights <= 0.0)) return("'weights' must all be positive")
           }
           if (!is.null(object@offset) && length(object@offset) != numObservations) return("'offset' must be null or have length equal to that of 'y'")
           if (!anyNA(object@n.cuts) && length(object@n.cuts) != ncol(object@x)) return("length of 'n.cuts' must equal number of columns in 'x'")
           
           if (!is.na(object@sigma) && object@sigma <= 0.0) return("'sigma' must be positive")
         })

## this shouldn't ever get created, used, modified, whathaveyou
setClass("dbartsState",
         slots =
         list(fit.tree    = "numeric",
              fit.total   = "numeric",
              fit.test    = "numeric",
              sigma       = "numeric",
              runningTime = "numeric",
              trees       = "character"))
              

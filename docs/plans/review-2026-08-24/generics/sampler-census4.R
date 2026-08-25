suppressPackageStartupMessages(library(dbarts))
grDevices::pdf(NULL)
OUT <- Sys.getenv("D")
csv <- file.path(OUT, "out2", "sampler-census4.csv")
esc <- function(s){s<-as.character(s);s[is.na(s)]<-"";s<-gsub('"','""',s,fixed=TRUE);s<-gsub("[\r\n]+"," ",s);paste0('"',s,'"')}
conn <- file(csv, "wt"); writeLines("sampler,phase,method,outcome,detail", conn)
raw <- c("^'arg' should be one of","unused argument","subscript out of bounds","could not find function",
 "is missing, with no default","missing value where TRUE/FALSE needed","non-numeric argument",
 "invalid subscript","object '.*' not found","attempt to apply non-function",
 "[$] operator is invalid for atomic vectors","invalid 'type'","NA/NaN/Inf",
 "unable to find an inherited method","no applicable method","incorrect number of dimensions",
 "non-conformable","cannot coerce","undefined columns selected","is not a function",
 "matched by multiple actual arguments","argument \"[a-zA-Z.]+\" is missing")
cls <- function(m){for(p in raw) if(grepl(p,m,perl=TRUE)) return("error-without-reason"); "refused"}
desc <- function(v){ if(is.null(v)) return("NULL"); d<-tryCatch(dim(v),error=function(e)NULL)
  paste(paste(class(v),collapse="/"), if(!is.null(d)) paste(d,collapse="x") else paste0("len",length(v))) }
SAMP <- ""; PHASE <- ""
cm <- function(method, expr) {
  o <- withCallingHandlers(tryCatch(list(ok=TRUE,v=force(expr),m=NA_character_),
        error=function(e) list(ok=FALSE,v=NULL,m=conditionMessage(e))),
        warning=function(w) invokeRestart("muffleWarning"), message=function(m) invokeRestart("muffleMessage"))
  writeLines(paste(esc(SAMP),esc(PHASE),esc(method),
    esc(if(o$ok) "accepted" else cls(o$m)), esc(if(o$ok) desc(o$v) else o$m), sep=","), conn); flush(conn)
  invisible(o)
}

set.seed(20260824); n <- 40L; nTest <- 10L
x <- matrix(runif(n*3L),n,3L); colnames(x) <- c("x1","x2","x3")
xTest <- matrix(runif(nTest*3L),nTest,3L); colnames(xTest) <- colnames(x)
yG <- as.numeric(2*x[,1]-x[,2]+rnorm(n,0,0.3))
counts <- t(rmultinom(n, size=10L, prob=c(0.4,0.35,0.25))); colnames(counts) <- c("a","b","c")
yBin <- rbinom(n,1L,plogis(2*(x[,1]-0.5)))
wInt <- rep(2L, n)
ctrl <- function(...) do.call(dbartsControl, utils::modifyList(
  list(n.chains=1L,n.threads=1L,n.trees=4L,n.burn=3L,n.samples=5L,keepTrees=TRUE,keepTrainingFits=TRUE), list(...)))

mk <- list(
  multinomial = function() dbarts(dbartsData(x, counts=counts), family="multinomial", control=ctrl(), verbose=FALSE),
  amplitude   = function() dbarts(x, yG, control=ctrl(),
                    forests=list(forest(n.trees=3L), forest(basis=x[,2,drop=FALSE], n.trees=3L)), verbose=FALSE),
  hetero      = function() dbarts(x, yG, control=ctrl(), variance=~x1, verbose=FALSE),
  wlogistic   = function() dbarts(x, yBin, family="logistic", weights=wInt, control=ctrl(), verbose=FALSE)
)

for (sname in names(mk)) {
  SAMP <- sname
  s <- tryCatch(mk[[sname]](), error=function(e) e)
  if (inherits(s,"error")) { PHASE<-"build"; cm("construct", stop(conditionMessage(s))); next }
  PHASE <- "build"; cm("construct", s)
  rr <- tryCatch(s$run(3L,5L), error=function(e) e)
  PHASE <- "run"; cm("run(3,5)", if (inherits(rr,"error")) stop(conditionMessage(rr)) else rr)
  ## ---- reads first: the tree store is untouched here ----
  PHASE <- "read"
  cm("predict",            s$predict(xTest))
  cm("predictForests",     s$predictForests(xTest))
  cm("getTrees",           s$getTrees())
  cm("printTrees",         s$printTrees())
  cm("plotTree",           s$plotTree(1L,1L,1L))
  cm("getDispersion",      s$getDispersion())
  cm("getFitsWithoutOffset", s$getFitsWithoutOffset())
  cm("getLatents",         s$getLatents(rr))
  cm("getSigmas",          s$getSigmas(rr))
  cm("getSumsOfSquaredResiduals", s$getSumsOfSquaredResiduals(rr))
  cm("getForestFits",      s$getForestFits(1L))
  cm("getForestAmplitudes", s$getForestAmplitudes())
  cm("getForestVariableCounts", s$getForestVariableCounts(1L))
  cm("getCalibration",     s$getCalibration(1L))
  cm("copy",               s$copy(shallow=TRUE))
  cm("show",               s$show())
  cm("extract(S3)",        extract(s))
  cm("plotTree(S3)",       plotTree(s, 1L, 1L, 1L))
  cm("storeState",         s$storeState())
  cm("getPointer",         s$getPointer())
  ## ---- mutators ----
  PHASE <- "mutate"
  cm("setControl",         s$setControl(s$control))
  cm("setModel",           s$setModel(s$model))
  cm("setData",            s$setData(s$data))
  cm("setResponse",        s$setResponse(rnorm(n)))
  cm("setOffset",          s$setOffset(rep(0,n)))
  cm("setWeights",         s$setWeights(runif(n,0.5,1.5)))
  cm("setCounts",          s$setCounts(counts))
  cm("setCategoryOffset",  s$setCategoryOffset(rep(0,n)))
  cm("setCategoryTestOffset", s$setCategoryTestOffset(rep(0,nTest)))
  cm("setActiveRows",      s$setActiveRows(rep(TRUE,n)))
  cm("setForestWeights",   s$setForestWeights(1L, runif(n)))
  cm("setForestBasis",     s$setForestBasis(1L, ~x1))
  cm("setSigma",           s$setSigma(1.0))
  cm("setPredictor",       s$setPredictor(x))
  cm("setCutPoints",       s$setCutPoints(sort(runif(20)), column=1L))
  cm("setTestPredictorAndOffset", s$setTestPredictorAndOffset(xTest, NULL))
  cm("setTestPredictor",   s$setTestPredictor(xTest[,1L], column=1L))
  cm("setTestOffset",      s$setTestOffset(rep(0,nTest)))
  cm("setCalibration",     s$setCalibration(prior.scale=1.0))
  cm("installTrees",       { s2 <- s$copy(shallow=TRUE); s$installTrees(s2) })
  cm("sampleTreesFromPrior", s$sampleTreesFromPrior())
  cm("sampleNodeParametersFromPrior", s$sampleNodeParametersFromPrior())
  cm("growFromRoot",       s$growFromRoot(1L))
  cm("setState",           s$setState(s$storeState()))
  cm("adoptPointer",       s$adoptPointer(s$getPointer()))
  cm("run(0,2) after mutation", s$run(0L,2L))
}
close(conn); cat("census4 written\n")

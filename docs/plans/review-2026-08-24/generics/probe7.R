suppressPackageStartupMessages(library(dbarts))
set.seed(20260824); n<-40L; nTest<-10L
x <- matrix(runif(n*3L),n,3L); colnames(x)<-c("x1","x2","x3")
xTest <- matrix(runif(nTest*3L),nTest,3L); colnames(xTest)<-colnames(x)
yG <- as.numeric(2*x[,1]-x[,2]+rnorm(n,0,0.3)); zA <- rbinom(n,1L,0.5)
ct <- dbartsControl(n.chains=1L,n.threads=1L,n.trees=4L,n.burn=3L,n.samples=5L,keepTrees=TRUE,keepTrainingFits=TRUE)
mk <- function() { s <- dbarts(x, yG, control=ct, forests=list(forest(n.trees=3L), forest(basis=cbind(1-zA,zA), n.trees=3L)), verbose=FALSE); invisible(s$run(3L,5L)); s }
p <- function(lab, ex){ r <- tryCatch(withCallingHandlers(eval(ex), warning=function(w) invokeRestart("muffleWarning")),
    error=function(e) structure(conditionMessage(e), class="err"))
  if (inherits(r,"err")) cat(sprintf("%-56s ERROR: %s\n", lab, substr(as.character(r),1,110)))
  else cat(sprintf("%-56s OK %s %s\n", lab, paste(class(r),collapse="/"), if(!is.null(dim(r))) paste(dim(r),collapse="x") else paste0("len",length(r)))) }
cat("=== isolate: setForestBasis then state round-trip ===\n")
s1 <- mk(); s1$setForestBasis(2L, cbind(1-zA, zA))          # width-preserving
p("width-preserving basis; setState(state)", quote(s1$setState(s1$state)))
s2 <- mk(); s2$setForestBasis(1L, x[,3,drop=FALSE])          # forest 1 gains a basis (width 1)
p("forest1 gains basis; setState(state)",    quote(s2$setState(s2$state)))
s3 <- mk(); s3$setForestBasis(2L, x[,3,drop=FALSE])          # forest 2 narrows 2 -> 1
p("forest2 narrows 2->1; setState(state)",   quote(s3$setState(s3$state)))
s4 <- mk(); s4$setForestBasis(2L, x[,3,drop=FALSE]); s4$storeState()
p("forest2 narrows, explicit storeState",    quote(s4$setState(s4$state)))
s5 <- mk(); s5$setForestWeights(1L, runif(n))
p("setForestWeights; setState(state)",       quote(s5$setState(s5$state)))
s6 <- mk(); s6$setPredictor(x)
p("setPredictor; setState(state)",           quote(s6$setState(s6$state)))
cat("\n=== is a narrowed-basis sampler serializable? ===\n")
s7 <- mk(); s7$setForestBasis(2L, x[,3,drop=FALSE]); s7$storeState()
saveRDS(s7, file.path(Sys.getenv("D"),"s7.rds")); cat("saved s7\n")
cat("\n=== hurdle component storeState recipe ===\n")
yH <- ifelse(runif(n)<0.3, 0, rlnorm(n))
fH <- bart2(x, yH, family="hurdle.lognormal", n.trees=4L,n.samples=6L,n.burn=4L,n.chains=1L,n.threads=1L,
            keepTrees=TRUE,keepTrainingFits=TRUE,verbose=FALSE,seed=1L)
cat("occupancy names:", paste(names(fH$occupancy),collapse=","), "\n")
fH$occupancy$fit$storeState(); fH$positive$fit$storeState()
saveRDS(fH, file.path(Sys.getenv("D"),"fH2.rds")); cat("saved fH2\n")

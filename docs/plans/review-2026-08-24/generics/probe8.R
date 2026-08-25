suppressPackageStartupMessages(library(dbarts))
set.seed(20260824); n<-40L
x <- matrix(runif(n*3L),n,3L); colnames(x)<-c("x1","x2","x3")
yG <- as.numeric(2*x[,1]-x[,2]+rnorm(n,0,0.3)); zA <- rbinom(n,1L,0.5)
ct <- dbartsControl(n.chains=1L,n.threads=1L,n.trees=4L,n.burn=3L,n.samples=5L,keepTrees=TRUE,keepTrainingFits=TRUE)
mk <- function() { s <- dbarts(x, yG, control=ct, forests=list(forest(n.trees=3L), forest(basis=cbind(1-zA,zA), n.trees=3L)), verbose=FALSE); invisible(s$run(3L,5L)); s }
sA <- mk(); sA$setForestBasis(2L, x[,3,drop=FALSE])
stA <- sA$state
sB <- mk(); sB$setForestBasis(2L, x[,3,drop=FALSE]); sB$storeState(); stB <- sB$state
cat("promise state class:", class(stA), " length:", length(stA), "\n")
cat("store   state class:", class(stB), " length:", length(stB), "\n")
show1 <- function(st, lab) { e <- st[[1L]]; cat(lab, ": names=", paste(names(e), collapse=","), "\n")
  for (nm in names(e)) { v <- e[[nm]]; cat("   ", nm, ":", paste(class(v),collapse="/"),
    if(!is.null(dim(v))) paste(dim(v),collapse="x") else paste0("len",length(v)), "\n") } }
show1(stA, "promise"); show1(stB, "storeState")
cat("identical:", identical(stA, stB), "\n")
cat("setState(promise):", tryCatch({sA$setState(stA); "OK"}, error=function(e) conditionMessage(e)), "\n")
cat("setState(store)  :", tryCatch({sB$setState(stB); "OK"}, error=function(e) conditionMessage(e)), "\n")
cat("cross: sA$setState(stB):", tryCatch({sA$setState(stB); "OK"}, error=function(e) conditionMessage(e)), "\n")
# was the promise forced before the basis change?
sC <- mk(); invisible(sC$state); sC$setForestBasis(2L, x[,3,drop=FALSE])
cat("forced-before-change setState:", tryCatch({sC$setState(sC$state); "OK"}, error=function(e) conditionMessage(e)), "\n")

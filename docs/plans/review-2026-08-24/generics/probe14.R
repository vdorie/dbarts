suppressPackageStartupMessages({library(dbarts); library(survival)})
set.seed(21); n<-60L; nTest<-6L
x <- matrix(runif(n*2L),n,2L); colnames(x)<-c("x1","x2")
xTest <- matrix(runif(nTest*2L),nTest,2L); colnames(xTest)<-colnames(x)
tt <- rexp(n, rate=exp(-0.4*(x[,1]-0.5))); st <- rbinom(n,1L,0.85)
yS <- Surv(tt, st)
p <- function(lab, ex){ r <- tryCatch(withCallingHandlers(eval(ex),warning=function(w) invokeRestart("muffleWarning")),
    error=function(e) structure(conditionMessage(e),class="err"))
  if (inherits(r,"err")) cat(sprintf("%-52s ERROR: %s\n", lab, substr(as.character(r),1,110)))
  else cat(sprintf("%-52s OK %s %s\n", lab, paste(class(r),collapse="/"), if(!is.null(dim(r))) paste(dim(r),collapse="x") else paste0("len",length(r)))) }
fH <- bart2(x, yS, family="hazard", n.trees=4L,n.samples=6L,n.burn=4L,n.chains=2L,n.threads=1L,
            keepTrees=TRUE,keepTrainingFits=TRUE,verbose=FALSE,seed=1L)
cat("hazard fit: periods len", length(fH$periods), " family:", fH$family, " x cols:", paste(colnames(fH$fit$data@x),collapse=","), "\n")
p("survivalProbabilities(fH)  [training, default times]", quote(survivalProbabilities(fH)))
p("survivalProbabilities(fH, times=c(0.5,1))",  quote(survivalProbabilities(fH, times=c(0.5,1))))
p("survivalProbabilities(fH, newdata=xTest)",   quote(survivalProbabilities(fH, newdata=xTest)))
p("survivalProbabilities(fH, times=c(.5,1), newdata=xTest)", quote(survivalProbabilities(fH, times=c(0.5,1), newdata=xTest)))
p("survivalProbabilities(fH, combineChains=FALSE)", quote(survivalProbabilities(fH, combineChains=FALSE)))
fA <- bart2(x, yS, family="aft", n.trees=4L,n.samples=6L,n.burn=4L,n.chains=2L,n.threads=1L,
            keepTrees=TRUE,keepTrainingFits=TRUE,verbose=FALSE,seed=1L)
p("survivalProbabilities(fA, times=c(.5,1))",   quote(survivalProbabilities(fA, times=c(0.5,1))))
p("survivalProbabilities(fA, times=c(.5,1), combineChains=FALSE)", quote(survivalProbabilities(fA, times=c(0.5,1), combineChains=FALSE)))
p("survivalProbabilities(fA, times=c(.5,1), newdata=xTest)", quote(survivalProbabilities(fA, times=c(0.5,1), newdata=xTest)))
p("survivalProbabilities(fA)  [times missing]", quote(survivalProbabilities(fA)))
grp <- factor(sample(letters[1:3],n,TRUE)); grpT <- factor(sample(letters[1:3],nTest,TRUE), levels=letters[1:3])
df <- data.frame(x1=x[,1], x2=x[,2])
fR <- rbart_vi(yS ~ x1 + x2, df, group.by=grp, n.trees=4L,n.samples=30L,n.burn=10L,n.chains=1L,n.threads=1L,
               keepTrees=TRUE, verbose=FALSE)
p("survivalProbabilities(fR, times=c(.5,1))",  quote(survivalProbabilities(fR, times=c(0.5,1))))
p("survivalProbabilities(fR, times, newdata, group.by)", quote(survivalProbabilities(fR, times=c(0.5,1), newdata=data.frame(x1=xTest[,1],x2=xTest[,2]), group.by=grpT)))
cat("\n=== aft fit call spelling in the keepTrees error ===\n")
fA0 <- bart2(x, yS, family="aft", n.trees=4L,n.samples=6L,n.burn=4L,n.chains=1L,n.threads=1L,keepTrees=FALSE,keepSampler=FALSE,verbose=FALSE,seed=1L)
cat("call:", deparse(fA0$call)[1], "\n")
p("predict(fA0, xTest)", quote(predict(fA0, xTest)))

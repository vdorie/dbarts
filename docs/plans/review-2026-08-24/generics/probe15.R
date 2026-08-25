suppressPackageStartupMessages({library(dbarts); library(survival)})
set.seed(21); n<-60L; nTest<-4L
xn <- matrix(runif(n*2L),n,2L)                       # unnamed
xc <- xn; colnames(xc) <- c("x1","x2")               # named
tt <- rexp(n, rate=exp(-0.4*(xn[,1]-0.5))); st <- rbinom(n,1L,0.85)
yS <- Surv(tt, st)
mk <- function(X) bart2(X, yS, family="hazard", n.trees=4L,n.samples=6L,n.burn=4L,n.chains=1L,
        n.threads=1L,keepTrees=TRUE,keepTrainingFits=TRUE,verbose=FALSE,seed=1L)
p <- function(lab, ex){ r <- tryCatch(withCallingHandlers(eval(ex),warning=function(w) invokeRestart("muffleWarning")),
    error=function(e) structure(conditionMessage(e),class="err"))
  if (inherits(r,"err")) cat(sprintf("%-56s ERROR: %s\n", lab, substr(as.character(r),1,110)))
  else cat(sprintf("%-56s OK %s\n", lab, paste(dim(r),collapse="x"))) }
fN <- mk(xn); fC <- mk(xc)
cat("unnamed-x fit periods:", length(fN$periods), " named-x fit periods:", length(fC$periods), "\n")
p("UNNAMED x: survivalProbabilities(fit)",      quote(survivalProbabilities(fN)))
p("UNNAMED x: survivalProbabilities(fit, newdata=matrix)", quote(survivalProbabilities(fN, newdata=xn[1:nTest,,drop=FALSE])))
p("NAMED   x: survivalProbabilities(fit)",      quote(survivalProbabilities(fC)))
p("NAMED   x: survivalProbabilities(fit, newdata=named matrix)", quote(survivalProbabilities(fC, newdata=xc[1:nTest,,drop=FALSE])))
p("NAMED   x: survivalProbabilities(fit, newdata=data.frame)", quote(survivalProbabilities(fC, newdata=as.data.frame(xc[1:nTest,,drop=FALSE]))))
# formula-interface hazard fit (always named)
df <- data.frame(x1=xn[,1], x2=xn[,2])
fF <- bart2(yS ~ x1 + x2, df, family="hazard", n.trees=4L,n.samples=6L,n.burn=4L,n.chains=1L,
       n.threads=1L,keepTrees=TRUE,keepTrainingFits=TRUE,verbose=FALSE,seed=1L)
p("FORMULA hazard: survivalProbabilities(fit)", quote(survivalProbabilities(fF)))
p("FORMULA hazard: survivalProbabilities(fit, newdata=df)", quote(survivalProbabilities(fF, newdata=df[1:nTest,])))

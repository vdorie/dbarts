suppressMessages(library(dbarts))
say <- function(id, ...) cat("\n### ", id, ": ", ..., "\n", sep = "")
tc <- function(expr) tryCatch(expr, error = function(e) paste0("ERROR: ", conditionMessage(e)))
set.seed(99); n <- 60
x <- matrix(runif(n*3), n, 3); colnames(x) <- c("x1","x2","x3")
y <- 2*x[,1] - x[,2] + rnorm(n, 0, .5); yM <- factor(sample(c("a","b","c"), n, TRUE))
ctl <- dbartsControl(n.trees=3L, n.samples=5L, n.burn=3L, n.chains=1L, verbose=FALSE)
bs <- cbind(x[,3], 1 - x[,3])

say("G10", "dbartsSampler $state staleness after a mutation")
s <- dbarts(x, y, control=ctl, forests=list(forest(), forest(basis=bs)), verbose=FALSE)
s$run()
g0 <- length(s$state[[1L]]$glue)
r2 <- tc(s$setForestBasis(2L, matrix(1, n, 1L)))
cat("  setForestBasis(2, 1-col) -> ", substr(if (is.character(r2)) r2 else "ok",1,70), "\n")
cat("  glue length before mutation = ", g0, " ; s$state read AFTER mutation = ",
    tc(length(s$state[[1L]]$glue)), "\n", sep="")
s$storeState()
cat("  glue length after an explicit storeState() = ", length(s$state[[1L]]$glue), "\n")
s2 <- dbarts(x, y, control=ctl, forests=list(forest(), forest(basis=bs)), verbose=FALSE)
s2$run(); st <- s2$state
tc(s2$setForestBasis(2L, matrix(1, n, 1L)))
cat("  s2$setState(<stale state>) -> ", substr(tc({ s2$setState(st); "ACCEPTED" }),1,100), "\n")
s3 <- dbarts(x, y, control=ctl, forests=list(forest(), forest(basis=bs)), verbose=FALSE)
s3$run(); tc(s3$setForestBasis(2L, matrix(1, n, 1L))); s3$storeState()
cat("  storeState() first, then setState -> ", substr(tc({ s3$setState(s3$state); "ACCEPTED" }),1,100), "\n")

say("M7", "dbartsSampler-class.Rd:308 - setState(storeState()) round trip")
s4 <- dbarts(x, y, control=ctl, verbose=FALSE); s4$run()
cat("  s4$setState(s4$storeState()) -> ", substr(tc({ s4$setState(s4$storeState()); "ACCEPTED" }),1,100), "\n")
cat("  storeState() returns: ", tc(paste(class(s4$storeState()), collapse=",")), "\n")

say("M8", "setForestBasis(k, ~var) formula scope")
s5 <- dbarts(x, y, control=ctl, forests=list(forest(), forest(basis=bs)), verbose=FALSE)
cat("  s5$setForestBasis(2L, ~x3) -> ", substr(tc({ s5$setForestBasis(2L, ~x3); "ACCEPTED" }),1,100), "\n")
cat("  dbarts(forests=list(forest(), forest(basis=~x3))) at top level -> ",
    substr(tc({ dbarts(x, y, control=ctl, forests=list(forest(), forest(basis=~x3)), verbose=FALSE); "ACCEPTED" }),1,100), "\n")

say("M4", "setForestBasis vs setForestWeights on an amplitude-free sampler")
s6 <- dbarts(x, y, control=ctl, verbose=FALSE)
cat("  setForestBasis(1L, m)   -> ", substr(tc({ s6$setForestBasis(1L, matrix(1,n,1L)); "ACCEPTED" }),1,100), "\n")
cat("  setForestWeights(rep(1,n)) -> ", substr(tc({ s6$setForestWeights(rep(1,n)); "ACCEPTED" }),1,100), "\n")

say("E8", "per-category offset spellings")
oc <- matrix(0.1, n, 3)
cat("  bart2(offset=<nxK>)           -> ", substr(tc({ bart2(x, yM, family="multinomial", offset=oc, n.trees=3L, n.samples=3L, n.burn=2L, verbose=FALSE); "ACCEPTED" }),1,95), "\n")
cat("  bart2(offset.category=<nxK>)  -> ", substr(tc({ bart2(x, yM, family="multinomial", offset.category=oc, n.trees=3L, n.samples=3L, n.burn=2L, verbose=FALSE); "ACCEPTED" }),1,95), "\n")
cat("  dbarts(offset.category=<nxK>) -> ", substr(tc({ dbarts(x, yM, family="multinomial", offset.category=oc, control=ctl, verbose=FALSE); "ACCEPTED" }),1,95), "\n")
cat("  dbarts(offset=<nxK>)          -> ", substr(tc({ dbarts(x, yM, family="multinomial", offset=oc, control=ctl, verbose=FALSE); "ACCEPTED" }),1,95), "\n")
cnt <- model.matrix(~ yM - 1); storage.mode(cnt) <- "integer"
cat("  dbartsData(counts=, offset.category=) -> ", substr(tc({ dbartsData(x, counts=cnt, offset.category=oc); "ACCEPTED" }),1,95), "\n")
cat("  dbarts(<that data>, multinomial) -> ", substr(tc({ dbarts(dbartsData(x, counts=cnt, offset.category=oc), family="multinomial", control=ctl, verbose=FALSE); "ACCEPTED" }),1,95), "\n")

say("E19", "defaultNodeScale vs defaultAmplitudePriorScale")
for (fam in c("gaussian","probit","logistic","hazard","student","multinomial","nbinom")) {
  a <- tc(dbarts:::defaultNodeScale(fam)); b <- tc(dbarts:::defaultAmplitudePriorScale(fam))
  cat(sprintf("  %-12s nodeScale=%-14s amplitude=%s\n", fam,
      if (is.null(a)) "NULL" else paste(a, collapse=","),
      substr(if (is.character(b)) b else paste(b, collapse=","),1,60)))
}

say("E14", "makeind(all=)")
df <- data.frame(a = factor(c("p","q","r","p")), b = 1:4)
cat("  identical(makeind(df, all=FALSE), makeind(df)) = ", identical(makeind(df, all=FALSE), makeind(df)), "\n")
cat("  identical(makeind(df, all=TRUE),  makeind(df)) = ", identical(makeind(df, all=TRUE), makeind(df)), "\n")
cat("  body(makeind): ", paste(deparse(body(makeind)), collapse=" "), "\n")

say("E13", "dbartsDrawLatents sigma guard")
fv <- bart2(x, as.numeric(y > median(y)), n.trees=3L, n.samples=5L, n.burn=3L, n.chains=1L, verbose=FALSE, seed=8L)
fitv <- as.numeric(fitted(fv, type="bart")); yb <- as.numeric(y > median(y))
cat("  probit, sigma omitted -> ", substr(tc({v<-dbartsDrawLatents("probit", fit=fitv, y=yb); paste("numeric len", length(v))}),1,95), "\n")
cat("  probit, sigma = 1     -> ", substr(tc({v<-dbartsDrawLatents("probit", fit=fitv, y=yb, sigma=1); paste("numeric len", length(v))}),1,95), "\n")
cat("  aft, sigma omitted    -> ", substr(tc({v<-dbartsDrawLatents("aft", fit=fitv, y=cbind(abs(y)+.1, rbinom(n,1,.7))); paste("numeric len", length(v))}),1,95), "\n")
cat("  gaussian (R gate)     -> ", substr(tc({dbartsDrawLatents("gaussian", fit=fitv, y=y); "ACCEPTED"}),1,95), "\n")

say("M2", "residuals(fit, sample=)")
fB <- bart2(x, y, n.trees=3L, n.samples=6L, n.burn=3L, n.chains=1L, keepTrees=TRUE, verbose=FALSE, seed=5L)
cat("  residuals(fB, sample='train') -> ", substr(tc({residuals(fB, sample="train"); "ACCEPTED"}),1,95), "\n")

say("M3", "plotTree partial match")
pdf(NULL)
cat("  plotTree(fB, treeNum=1L, sample=2L) -> ", substr(tc({plotTree(fB, treeNum=1L, sample=2L); "DREW"}),1,95), "\n")
dev.off()

say("M1", "sample= validation message, bart vs own-class")
fO <- bart2(x, ordered(sample(1:3,n,TRUE)), family="ordinal", n.trees=3L, n.samples=6L, n.burn=3L, n.chains=1L, keepTrees=TRUE, verbose=FALSE, seed=5L)
cat("  extract(fB, sample='zzz') -> ", substr(tc(extract(fB, sample="zzz")),1,95), "\n")
cat("  extract(fO, sample='zzz') -> ", substr(tc(extract(fO, sample="zzz")),1,95), "\n")

say("M10", "names(fit) lists NULL components")
cat("  names(fB): ", paste(names(fB), collapse=","), "\n")
cat("  'yhat.test' %in% names(fB) = ", "yhat.test" %in% names(fB), " ; is.null(fB$yhat.test) = ", is.null(fB$yhat.test), "\n")

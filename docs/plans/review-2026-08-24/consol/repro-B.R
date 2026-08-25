suppressMessages(library(dbarts)); suppressMessages(library(survival))
say <- function(id, ...) cat("\n### ", id, ": ", ..., "\n", sep = "")
tc <- function(expr) tryCatch(expr, error = function(e) paste0("ERROR: ", conditionMessage(e)))
dm <- function(x) if (is.character(x)) x else paste(class(x)[1], paste(dim(x), collapse="x"))
set.seed(99)
n <- 60; x <- matrix(runif(n*3), n, 3); colnames(x) <- c("x1","x2","x3")
y <- 2*x[,1] - x[,2] + rnorm(n, 0, .5); xt <- matrix(runif(8*3), 8, 3); colnames(xt) <- colnames(x)
yb <- as.numeric(y > median(y)); yM <- factor(sample(c("a","b","c"), n, TRUE))
yO <- ordered(sample(1:3, n, TRUE)); yC <- rpois(n, 3)
yH <- ifelse(runif(n) < .5, 0, exp(rnorm(n)))
mk <- function(...) bart2(x, ..., n.trees=3L, n.samples=6L, n.burn=3L, n.chains=2L, keepTrees=TRUE, verbose=FALSE, seed=5L)
fB <- mk(y); fM <- mk(yM, family="multinomial"); fO <- mk(yO, family="ordinal")
fN <- mk(yC, family="nbinom"); fH <- mk(yH, family="hurdle.lognormal")
cat("classes:", class(fB)[1], class(fM)[1], class(fO)[1], class(fN)[1], class(fH)[1], "\n")

## G3/G4 argument swallow
say("G3/G4", "silently swallowed arguments on the own-class generics")
chk <- function(lbl, a, b) cat(sprintf("  %-46s identical-to-without = %s\n", lbl,
  tryCatch(identical(a, b), error=function(e) paste("ERR", conditionMessage(e)))))
chk("extract(fO, combineChains=FALSE) vs plain", tc(extract(fO, combineChains=FALSE)), tc(extract(fO)))
chk("extract(fM, combineChains=FALSE) vs plain", tc(extract(fM, combineChains=FALSE)), tc(extract(fM)))
chk("extract(fN, combineChains=FALSE) vs plain", tc(extract(fN, combineChains=FALSE)), tc(extract(fN)))
cat("  dims: extract(fO) ", dm(tc(extract(fO))), " ; predict(fO,xt,combineChains=FALSE) ",
    dm(tc(predict(fO, xt, combineChains=FALSE))), "\n")
chk("fitted(fM, sample='test') vs plain",   tc(fitted(fM, sample="test")), tc(fitted(fM)))
chk("fitted(fO, ci.level=0.9) vs plain",    tc(fitted(fO, ci.level=0.9)),  tc(fitted(fO)))
chk("predict(fN, xt, ci.level=0.9) vs plain", tc(predict(fN, xt, ci.level=0.9)), tc(predict(fN, xt)))
chk("residuals(fM, type='bart') vs plain",  tc(residuals(fM, type="bart")), tc(residuals(fM)))
chk("extract(fM, forest=1L) vs plain",      tc(extract(fM, forest=1L)),    tc(extract(fM)))
chk("summary(fM, vars='sigma') vs plain",   tc(summary(fM, vars="sigma")), tc(summary(fM)))
cat("  contrast bart: extract(fB, combineChains=FALSE) ", dm(tc(extract(fB, combineChains=FALSE))),
    " vs plain ", dm(tc(extract(fB))), "\n")
cat("  contrast bart: fitted(fB, ci.level=0.9) ", dm(tc(fitted(fB, ci.level=0.9))), "\n")
cat("  formals(extract.bartOrdinal): ", paste(names(formals(dbarts:::extract.bartOrdinal)), collapse=","), "\n")
cat("  formals(extract.bart):        ", paste(names(formals(dbarts:::extract.bart)), collapse=","), "\n")

## G5 n.chains
say("G5", "$n.chains presence by keepTrees/keepSampler")
for (kt in c(TRUE, FALSE)) for (ks in c(TRUE, FALSE)) {
  f <- bart2(x, y, n.trees=3L, n.samples=6L, n.burn=3L, n.chains=2L, keepTrees=kt,
             keepSampler=ks, verbose=FALSE, seed=6L)
  cat(sprintf("  keepTrees=%-5s keepSampler=%-5s -> n.chains %s ; has $fit %s\n", kt, ks,
      if (is.null(f$n.chains)) "ABSENT" else f$n.chains, !is.null(f$fit)))
}
cat("  siblings under keepTrees=TRUE: multinomial", !is.null(fM$n.chains),
    " ordinal", !is.null(fO$n.chains), " nbinom", !is.null(fN$n.chains),
    " hurdle", !is.null(fH$n.chains), " bart", !is.null(fB$n.chains), "\n")

## G6 extract trees on keepSampler-only
say("G6", "extract(type='trees') with keepTrees=FALSE, keepSampler=TRUE")
fks <- bart2(x, y, n.trees=3L, n.samples=6L, n.burn=3L, n.chains=1L, keepTrees=FALSE,
             keepSampler=TRUE, verbose=FALSE, seed=7L)
tks <- tc(extract(fks, type="trees"))
cat("  -> "); cat(if (is.data.frame(tks)) sprintf("data.frame %dx%d cols=%s\n", nrow(tks), ncol(tks),
    paste(names(tks), collapse=",")) else paste(tks, "\n"))
fkt <- bart2(x, y, n.trees=3L, n.samples=6L, n.burn=3L, n.chains=1L, keepTrees=TRUE, verbose=FALSE, seed=7L)
tkt <- extract(fkt, type="trees")
cat("  keepTrees=TRUE cols = ", paste(names(tkt), collapse=","), "\n")

## M5: chain column on a single-chain fit
say("M5", "trees table columns, 1 chain vs 2 chains")
f2c <- bart2(x, y, n.trees=3L, n.samples=6L, n.burn=3L, n.chains=2L, keepTrees=TRUE, verbose=FALSE, seed=7L)
cat("  1 chain: ", paste(names(tkt), collapse=","), "\n  2 chains: ",
    paste(names(extract(f2c, type="trees")), collapse=","), "\n")

## G7 plot
say("G7", "plot() on the own-class fits")
for (nm in c("fB","fM","fO","fN","fH")) {
  f <- get(nm); pdf(NULL)
  r <- tc({ plot(f); "DREW" }); dev.off()
  cat(sprintf("  %-3s (%-16s) -> %s\n", nm, class(f)[1], substr(r,1,80)))
  cat(sprintf("       getS3method('plot','%s') = %s\n", class(f)[1],
      !is.null(getS3method("plot", class(f)[1], optional=TRUE))))
}

## G8 loglik
say("G8", "extract(type='loglik') per class")
for (nm in c("fB","fM","fO","fN","fH")) {
  f <- get(nm); cat(sprintf("  %-16s -> %s\n", class(f)[1], substr(dm(tc(extract(f, type="loglik"))),1,90)))
}

## G9 hurdle reload
say("G9", "saveRDS/readRDS then predict, per class")
td <- tempfile(); dir.create(td)
for (nm in c("fB","fM","fO","fN","fH")) {
  f <- get(nm); pth <- file.path(td, paste0(nm, ".rds"))
  tc(if (!is.null(f$fit)) f$fit$storeState())
  saveRDS(f, pth)
  out <- system2("Rscript", c("-e", shQuote(sprintf(
    'suppressMessages(library(dbarts)); f <- readRDS("%s"); xt <- matrix(runif(24),8,3); colnames(xt) <- c("x1","x2","x3"); cat(tryCatch(paste(class(predict(f, xt))[1], paste(dim(predict(f, xt)), collapse="x")), error=function(e) paste("ERROR:", conditionMessage(e))))',
    pth))), stdout=TRUE, stderr=TRUE)
  cat(sprintf("  %-16s -> %s\n", class(f)[1], substr(paste(out, collapse=" "),1,110)))
}
cat("  hurdle names(): ", paste(names(fH), collapse=","), "\n")

## G10 stale $state
say("G10", "dbartsSampler $state after a mutation (Rd:328 'reading it forces the current state')")
s <- dbarts(x, y, control=dbartsControl(n.trees=3L, n.samples=5L, n.burn=3L, n.chains=1L, verbose=FALSE),
            forests=list(forest(), forest(basis=~x3)), verbose=FALSE)
r <- tc(s$run())
cat("  run -> ", substr(if (is.character(r)) r else "ok", 1, 70), "\n")
if (!is.character(r)) {
  g0 <- tc(length(s$state[[1L]]$glue))
  r2 <- tc(s$setForestBasis(2L, matrix(1, n, 1L)))
  cat("  setForestBasis(2, 1-col) -> ", substr(if (is.character(r2)) r2 else "ok",1,70), "\n")
  cat("  length(s$state[[1]]$glue) after mutation = ", tc(length(s$state[[1L]]$glue)),
      " (before = ", g0, ")\n", sep="")
  cat("  fresh storeState glue length = ", tc({ s$storeState(); length(s$state[[1L]]$glue) }), "\n")
}
s2 <- dbarts(x, y, control=dbartsControl(n.trees=3L, n.samples=5L, n.burn=3L, n.chains=1L, verbose=FALSE),
             forests=list(forest(), forest(basis=~x3)), verbose=FALSE)
s2$run(); st <- s2$state
tc(s2$setForestBasis(2L, matrix(1, n, 1L)))
cat("  s2$setState(stale state) -> ", substr(tc({ s2$setState(st); "ACCEPTED" }), 1, 90), "\n")

## E8 offset.category
say("E8", "per-category offset spellings")
oc <- matrix(0.1, n, 3)
cat("  bart2(offset=<nxK>)            -> ", substr(tc({ bart2(x, yM, family="multinomial", offset=oc, n.trees=3L, n.samples=3L, n.burn=2L, verbose=FALSE); "ACCEPTED" }),1,90), "\n")
cat("  bart2(offset.category=<nxK>)   -> ", substr(tc({ bart2(x, yM, family="multinomial", offset.category=oc, n.trees=3L, n.samples=3L, n.burn=2L, verbose=FALSE); "ACCEPTED" }),1,90), "\n")
cat("  dbarts(offset.category=<nxK>)  -> ", substr(tc({ dbarts(x, yM, family="multinomial", offset.category=oc, verbose=FALSE); "ACCEPTED" }),1,90), "\n")
cat("  dbarts(offset=<nxK>)           -> ", substr(tc({ dbarts(x, yM, family="multinomial", offset=oc, verbose=FALSE); "ACCEPTED" }),1,90), "\n")
cnt <- model.matrix(~ yM - 1); storage.mode(cnt) <- "integer"
cat("  dbartsData(counts=, offset.category=) -> ", substr(tc({ dbartsData(x, counts=cnt, offset.category=oc); "ACCEPTED" }),1,90), "\n")

## E19 defaultNodeScale
say("E19", "defaultNodeScale vs defaultAmplitudePriorScale on an unknown family")
for (fam in c("gaussian","probit","hazard","student","multinomial")) {
  a <- tc(dbarts:::defaultNodeScale(fam)); b <- tc(dbarts:::defaultAmplitudePriorScale(fam))
  cat(sprintf("  %-12s nodeScale=%-28s amplitude=%s\n", fam,
      if (is.null(a)) "NULL" else paste(a, collapse=","), substr(if (is.character(b)) b else paste(b, collapse=","),1,50)))
}

## E14 makeind(all=)
say("E14", "makeind(all=)")
df <- data.frame(a = factor(c("p","q","r","p")), b = 1:4)
cat("  identical(makeind(df, all=FALSE), makeind(df)) = ", identical(makeind(df, all=FALSE), makeind(df)), "\n")
cat("  identical(makeind(df, all=TRUE),  makeind(df)) = ", identical(makeind(df, all=TRUE), makeind(df)), "\n")

## corrections
say("CORRECTIONS", "E13, E4-rbart, G1 data.frame newdata")
fv <- bart2(x, yb, n.trees=3L, n.samples=5L, n.burn=3L, n.chains=1L, verbose=FALSE, seed=8L)
fitv <- as.numeric(fitted(fv, type="bart"))
cat("  drawLatents(probit, sigma omitted) -> ", substr(tc({v<-dbartsDrawLatents("probit", fit=fitv, y=yb); paste("numeric len", length(v))}),1,90), "\n")
cat("  drawLatents(probit, sigma = 1)     -> ", substr(tc({v<-dbartsDrawLatents("probit", fit=fitv, y=yb, sigma=1); paste("numeric len", length(v))}),1,90), "\n")
cat("  drawLatents(aft, sigma omitted)    -> ", substr(tc({v<-dbartsDrawLatents("aft", fit=fitv, y=abs(y)+.1); paste("numeric len", length(v))}),1,90), "\n")
g <- factor(sample(1:4, n, TRUE)); tt <- rexp(n, .3); ss <- rbinom(n, 1, .7)
cat("  rbart_vi(x, Surv, aft) [n.samples 40] -> ", substr(tc(rbart_vi(x, Surv(tt,ss), group.by=g, family="aft", n.samples=40L, n.burn=10L, n.trees=3L, n.chains=1L, verbose=FALSE)),1,110), "\n")
fhN <- bart2(x, Surv(tt,ss), family="hazard", n.trees=3L, n.samples=5L, n.burn=3L, n.chains=1L, keepTrees=TRUE, verbose=FALSE, seed=2L)
cat("  survivalProbabilities(named fit, newdata=<data.frame>) -> ",
    substr(dm(tc(survivalProbabilities(fhN, newdata=as.data.frame(xt)))),1,90), "\n")
cat("  survivalProbabilities(named fit, newdata=<named matrix>) -> ",
    substr(dm(tc(survivalProbabilities(fhN, newdata=xt))),1,90), "\n")

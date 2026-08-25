## Consolidator independent reproduction, part A: R surface claims.
suppressMessages(library(dbarts)); suppressMessages(library(survival))
cat("dbarts", as.character(packageVersion("dbarts")), "\n")
say <- function(id, ...) cat("\n### ", id, ": ", ..., "\n", sep = "")
tc <- function(expr) tryCatch(expr, error = function(e) paste0("ERROR: ", conditionMessage(e)),
                              warning = function(w) paste0("WARNING: ", conditionMessage(w)))
p <- function(x) { if (is.character(x) && length(x) == 1L) cat("  ", x, "\n") else
                   cat("  ", paste(utils::capture.output(str(x, max.level = 1L)), collapse = " | "), "\n") }

set.seed(99)
n <- 60; p2 <- 3
x <- matrix(runif(n * p2), n, p2); colnames(x) <- c("x1","x2","x3")
y <- 2 * x[,1] - x[,2] + rnorm(n, 0, 0.5)
xt <- matrix(runif(12 * p2), 12, p2); colnames(xt) <- colnames(x)
yb <- as.numeric(y > median(y))
yM <- factor(sample(c("a","b","c"), n, TRUE))
yO <- ordered(sample(1:3, n, TRUE))
yC <- rpois(n, 3)

## ---------- E1 / G2 / F3 : extract(type="trees", sample=) ----------
say("E1/G2/F3", "extract(type='trees', sample=)")
f <- bart2(x, y, n.trees = 3L, n.samples = 5L, n.burn = 3L, n.chains = 1L,
           keepTrees = TRUE, verbose = FALSE, seed = 1L)
tr0 <- tc(extract(f, type = "trees"))
cat("  omitted        -> "); p(if (is.data.frame(tr0)) sprintf("data.frame %dx%d cols=%s", nrow(tr0), ncol(tr0), paste(names(tr0), collapse=",")) else tr0)
tr1 <- tc(extract(f, type = "trees", sample = 1L))
cat("  sample=1L      -> "); p(if (is.data.frame(tr1)) sprintf("data.frame %dx%d", nrow(tr1), ncol(tr1)) else tr1)
cat("  sample='train' -> "); p(tc(suppressWarnings(extract(f, type = "trees", sample = "train"))))
cat("  positional     -> "); p(tc(suppressWarnings(extract(f, "trees", "train"))))
cat("  combineChains=FALSE -> "); p(tc(extract(f, type = "trees", combineChains = FALSE)))
cat("  forest=1L      -> "); p(tc(extract(f, type = "trees", forest = 1L)))
cat("  contribution=TRUE -> "); p(tc(extract(f, type = "trees", contribution = TRUE)))
cat("  getTrees formals: "); cat(paste(names(formals(f$fit$getTrees)), collapse=","), "\n")

## ---------- G1 : survivalProbabilities on named-column hazard fit ----------
say("G1", "survivalProbabilities on a hazard fit with named x columns")
set.seed(7)
tt <- rexp(n, 0.3); ss <- rbinom(n, 1, 0.7)
fhN <- tc(bart2(x, Surv(tt, ss), family = "hazard", n.trees = 3L, n.samples = 5L,
                n.burn = 3L, n.chains = 1L, keepTrees = TRUE, verbose = FALSE, seed = 2L))
cat("  named-matrix fit: "); p(if (is.character(fhN)) fhN else class(fhN))
if (!is.character(fhN)) { cat("  survivalProbabilities(named) -> "); p(tc(survivalProbabilities(fhN))) }
xu <- x; colnames(xu) <- NULL
fhU <- tc(bart2(xu, Surv(tt, ss), family = "hazard", n.trees = 3L, n.samples = 5L,
                n.burn = 3L, n.chains = 1L, keepTrees = TRUE, verbose = FALSE, seed = 2L))
if (!is.character(fhU)) { r <- tc(survivalProbabilities(fhU))
  cat("  survivalProbabilities(unnamed) -> "); p(if (is.character(r)) r else paste("array", paste(dim(r), collapse="x"))) }
xdf <- as.data.frame(x)
fhD <- tc(bart2(xdf, Surv(tt, ss), family = "hazard", n.trees = 3L, n.samples = 5L,
                n.burn = 3L, n.chains = 1L, keepTrees = TRUE, verbose = FALSE, seed = 2L))
if (!is.character(fhD)) { r <- tc(survivalProbabilities(fhD))
  cat("  survivalProbabilities(data.frame) -> "); p(if (is.character(r)) r else paste("array", paste(dim(r), collapse="x"))) }
cat("  tinytest hazard fixture has colnames? ")
cat(any(grepl("colnames", readLines(system.file("tinytest","test-hazard.R", package="dbarts"))[1:40])), "\n")

## ---------- E2 / E3 / F9 : bart() family redirects ----------
say("E2/E3/F9", "bart() family token refusals")
for (fam in c("gaussian","probit","logistic","aft","multinomial","ordinal","nbinom",
              "hurdle.lognormal","twopart","hazard","hazard.probit","hazard.logistic")) {
  m <- tc(bart(x, y, family = fam, ndpost = 3L, nskip = 2L, ntree = 3L, verbose = FALSE))
  cat(sprintf("  %-18s -> %s\n", fam, if (is.character(m)) substr(m,1,90) else "ACCEPTED"))
}
cat("  bartOwnClassFamilies = "); cat(paste(dbarts:::bartOwnClassFamilies, collapse=", "), "\n")

## ---------- E4 / F2 : dbartsData multi-column response row count ----------
say("E4/F2", "dbartsData with a multi-column response")
sv <- Surv(tt, ss)
cat("  NROW(x)=", NROW(x), " NROW(sv)=", NROW(sv), " dim(sv)=", paste(dim(sv), collapse="x"), "\n")
cat("  dbartsData(x, sv)            -> "); p(tc(dbartsData(x, sv)))
cat("  dbartsData(x, cbind(y,y))    -> "); p(tc(dbartsData(x, cbind(y, y))))
cat("  rbart_vi(x, sv, aft)         -> ");
g <- factor(sample(1:4, n, TRUE))
p(tc(rbart_vi(x, sv, group.by = g, family = "aft", n.samples = 3L, n.burn = 2L, n.trees = 3L, n.chains = 1L, verbose = FALSE)))
cat("  xbart(x, sv)                 -> "); p(tc(xbart(x, sv, n.reps = 1L, n.trees = 3L, n.samples = 3L, n.burn = c(2L,1L), verbose = FALSE)))
cat("  dbarts(x, sv, aft)           -> "); p(tc(class(dbarts(x, sv, family = "aft", control = dbartsControl(n.trees=3L, n.samples=3L, n.burn=2L, verbose=FALSE), verbose = FALSE))))

## ---------- E5 : combinechains shape ----------
say("E5", "combineChains result shape (bart.Rd:153 says nchain x ndpost x nobs when TRUE)")
f2t <- bart2(x, y, n.trees = 3L, n.samples = 6L, n.burn = 3L, n.chains = 2L,
             combineChains = TRUE, verbose = FALSE, seed = 3L)
f2f <- bart2(x, y, n.trees = 3L, n.samples = 6L, n.burn = 3L, n.chains = 2L,
             combineChains = FALSE, verbose = FALSE, seed = 3L)
cat("  combineChains=TRUE  yhat.train dim = ", paste(dim(f2t$yhat.train), collapse="x"), "\n")
cat("  combineChains=FALSE yhat.train dim = ", paste(dim(f2f$yhat.train), collapse="x"), "\n")
b2t <- bart(x, y, ndpost = 6L, nskip = 3L, ntree = 3L, nchain = 2L, combinechains = TRUE, verbose = FALSE)
b2f <- bart(x, y, ndpost = 6L, nskip = 3L, ntree = 3L, nchain = 2L, combinechains = FALSE, verbose = FALSE)
cat("  bart() TRUE  dim = ", paste(dim(b2t$yhat.train), collapse="x"), "\n")
cat("  bart() FALSE dim = ", paste(dim(b2f$yhat.train), collapse="x"), "\n")

## ---------- E6 : proposalprobs defaults ----------
say("E6", "proposal.probs live default (bart.Rd:165 lists birth_death, change, swap = .5,.1,.4)")
cat("  dbarts formals: "); print(eval(formals(dbarts)$proposal.probs))
cat("  bart2 formals : "); print(eval(formals(bart2)$proposal.probs))
cat("  bart formals  : "); print(eval(formals(bart)$proposalprobs))

## ---------- E7 / F1 : dbartsSpec vs dbarts on multinomial ----------
say("E7/F1", "dbartsSpec vs dbarts, family='multinomial' on a factor response")
ctl <- dbartsControl(n.trees = 3L, n.samples = 5L, n.burn = 3L, n.chains = 1L, verbose = FALSE)
dd <- tc(dbartsData(x, yM))
cat("  dbartsData(x, yM)                        -> "); p(if (is.character(dd)) dd else class(dd))
if (!is.character(dd)) { cat("  dbartsSpec(dd, ctl, 'multinomial')        -> "); p(tc(class(dbartsSpec(data = dd, control = ctl, family = "multinomial")))) }
cat("  dbarts(x, yM, family='multinomial')      -> ")
p(tc({ s <- dbarts(x, yM, family = "multinomial", control = ctl, verbose = FALSE); paste("ACCEPTED, family =", s$model@family) }))
cnt <- model.matrix(~ yM - 1); storage.mode(cnt) <- "integer"
ddc <- tc(dbartsData(x, counts = cnt))
if (!is.character(ddc)) { cat("  dbartsSpec(counts-built dd, 'multinomial')-> "); p(tc(class(dbartsSpec(data = ddc, control = ctl, family = "multinomial")))) }
cat("  dbartsSpec family vocabulary: "); cat(paste(eval(formals(dbartsSpec)$family), collapse=", "), "\n")

## ---------- E9 / E10 / E11 : cross-entry disagreements ----------
say("E9", "bad family AND a Surv response: which is reported first")
for (fn in c("dbarts","bart2","bart","rbart_vi","xbart")) {
  m <- switch(fn,
    dbarts   = tc(dbarts(x, sv, family = "zzz", control = ctl, verbose = FALSE)),
    bart2    = tc(bart2(x, sv, family = "zzz", n.trees=3L, n.samples=3L, n.burn=2L, verbose=FALSE)),
    bart     = tc(bart(x, sv, family = "zzz", ntree=3L, ndpost=3L, nskip=2L, verbose=FALSE)),
    rbart_vi = tc(rbart_vi(x, sv, group.by = g, family = "zzz", n.trees=3L, n.samples=3L, n.burn=2L, verbose=FALSE)),
    xbart    = tc(xbart(x, sv, family = "zzz", n.reps=1L, n.trees=3L, n.samples=3L, n.burn=c(2L,1L))))
  cat(sprintf("  %-9s -> %s\n", fn, substr(if (is.character(m)) m else "ACCEPTED", 1, 90)))
}
say("E10", "n.threads = c(1, 2)")
for (fn in c("dbarts","bart2","bart","rbart_vi","xbart")) {
  m <- switch(fn,
    dbarts   = tc(dbarts(x, y, control = dbartsControl(n.threads = c(1L,2L)), verbose = FALSE)),
    bart2    = tc(bart2(x, y, n.threads = c(1L,2L), n.trees=3L, n.samples=3L, n.burn=2L, verbose=FALSE)),
    bart     = tc(bart(x, y, nthread = c(1L,2L), ntree=3L, ndpost=3L, nskip=2L, verbose=FALSE)),
    rbart_vi = tc(rbart_vi(x, y, group.by = g, n.threads = c(1L,2L), n.trees=3L, n.samples=3L, n.burn=2L, verbose=FALSE)),
    xbart    = tc(xbart(x, y, n.threads = c(1L,2L), n.reps=1L, n.trees=3L, n.samples=3L, n.burn=c(2L,1L))))
  cat(sprintf("  %-9s -> %s\n", fn, substr(if (is.character(m)) m else "ACCEPTED", 1, 90)))
}
say("E11", "n.samples = 0")
cat("  dbartsControl(n.samples=0) -> "); p(tc(class(dbartsControl(n.samples = 0L))))
cat("  dbarts(n.samples=0)        -> "); p(tc(class(dbarts(x, y, control = dbartsControl(n.samples = 0L, n.trees=3L), verbose = FALSE))))
cat("  bart2(n.samples=0)         -> "); p(tc(bart2(x, y, n.samples = 0L, n.trees=3L, n.burn=2L, verbose=FALSE)))
cat("  xbart(n.samples=0)         -> "); p(tc(xbart(x, y, n.samples = 0L, n.reps=1L, n.trees=3L, n.burn=c(2L,1L))))
cat("  rbart_vi(n.samples=0)      -> "); p(tc(rbart_vi(x, y, group.by=g, n.samples = 0L, n.trees=3L, n.burn=2L, verbose=FALSE)))

## ---------- E12 : monotone spellings ----------
say("E12", "monotone token vocabulary")
for (s in list("+","increasing",1,"-","decreasing",-1,"inc","dec","0","INC","Increasing","up")) {
  m <- tc(dbarts(x, y, monotone = setNames(list(s), "x1"),
                 control = dbartsControl(n.trees=3L, n.samples=3L, n.burn=2L, verbose=FALSE), verbose = FALSE))
  cat(sprintf("  %-12s -> %s\n", paste0(deparse(s), collapse=""), substr(if (is.character(m)) m else "ACCEPTED", 1, 80)))
}

## ---------- E13 : dbartsDrawLatents sigma default ----------
say("E13", "dbartsDrawLatents(family='probit', sigma = <the formal default>)")
cat("  augmentation formals: "); print(formals(dbartsDrawLatents))
fb <- bart2(x, yb, n.trees = 3L, n.samples = 5L, n.burn = 3L, n.chains = 1L, keepTrees = TRUE, verbose = FALSE, seed = 4L)
cat("  omitted   -> "); p(tc({ v <- dbartsDrawLatents("probit", fit = fb, y = yb); paste("numeric len", length(v)) }))
cat("  sigma = 1 -> "); p(tc({ v <- dbartsDrawLatents("probit", fit = fb, y = yb, sigma = 1); paste("numeric len", length(v)) }))

suppressMessages(library(dbarts))
bcCal <- getFromNamespace("bartcoreForestCalibration", "dbarts")
bcSet <- getFromNamespace("bartcoreSetForestPriorScale", "dbarts")
bcRun <- getFromNamespace("bartcoreRun", "dbarts")
set.seed(11L); n <- 160L; p <- 3L
x <- matrix(runif(n * p), n, p)
ctl <- dbartsControl(n.trees = 50L, n.chains = 1L, n.threads = 1L,
                     updateState = FALSE, verbose = FALSE, seed = 271L)
keep <- c("prior.scale","prior.sd","prior.mean","response.scale")
calR <- function(s) s$getCalibration(1L)[1L, keep]
calC <- function(s) bcCal(s, 0L)[1L, keep]
yA <- seq(-2.5, 2.5, length.out = n); yB <- 4 * yA + 7
S <- 2.5
mkAft <- function(y, status, np = NULL) {
  args <- list(x, y, control = ctl)
  if (!is.null(np)) args$node.prior <- np
  s <- do.call(dbarts, args)
  c2 <- s$control; attr(c2, "bartcore.survival") <- as.numeric(status); s$control <- c2
  dbarts:::bartcoreSampler(s, family = "aft")
}
st <- rep(c(1, 1, 1, 0), length.out = n)
cat("--- aft (internal bartcore surface), default vs named ---\n")
print(rbind(`aft anchor  (default)` = calC(mkAft(yA, st)),
            `aft rebuild (default)` = calC(mkAft(yB, st)),
            `aft anchor  (named)`   = calC(mkAft(yA, st, dbartsPriors$normal(2, scale = S))),
            `aft rebuild (named)`   = calC(mkAft(yB, st, dbartsPriors$normal(2, scale = S)))))
mkGrp <- function(y, np = NULL) {
  c2 <- ctl
  attr(c2, "bartcore.groups") <- list(indices = rep_len(1:8, n), n.groups = 8L,
                                      prior = "gamma", rel.scale = 0.2, n.steps = 5L)
  args <- list(x, y, control = c2)
  if (!is.null(np)) args$node.prior <- np
  do.call(dbarts, args)
}
cat("\n--- grouped gaussian (R5 surface), default vs named ---\n")
print(rbind(`grp anchor  (default)` = calR(mkGrp(yA)),
            `grp rebuild (default)` = calR(mkGrp(yB)),
            `grp anchor  (named)`   = calR(mkGrp(yA, dbartsPriors$normal(2, scale = S))),
            `grp rebuild (named)`   = calR(mkGrp(yB, dbartsPriors$normal(2, scale = S)))))
cat("\n--- mid-chain setCalibration ---\n")
g <- mkGrp(yB); g$setCalibration(prior.scale = S)
cat("grouped after setCalibration: "); print(calR(g))
a <- mkAft(yB, st, NULL)
r <- tryCatch({ bcSet(a, 0L, S); calC(a) }, error = function(e) paste("REFUSED:", conditionMessage(e)))
cat("aft after setCalibration:     "); print(r)
cat("\ngrouped run ok: ", !is.null(g$run(20L, 5L)$sigma), "\n")

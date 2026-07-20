# forest-ranef-collapse-proto.R
#
# Isolation prototype for docs/design/forest-ranef-interweaving.md. It measures
# the mixing gain of the PRIMARY remedy (a joint / collapsed (f, b) move that
# marginalizes the group intercepts out of the mean-fit update) against the
# centered blocked Gibbs the engine runs today, on a TRACTABLE surrogate:
#
#   y_i = beta0 + beta1 x_i + b_{g(i)} + eps_i,   b_j ~ N(0, tau^2),
#         eps ~ N(0, sigma^2),   sigma known.
#
# The linear mean beta1 x plays the role of the BART forest f; x is clustered by
# group (x_i = group_center + N(0, s)) so beta1 x competes with b_g to explain
# each group's mean -- the identifiability confounding grouped-mixing.R shows
# dominates tau's mixing at small K. This is a caricature (f linear, so the
# collapse is closed form; the real forest collapse is far harder -- see the
# design note), but it isolates the SAME beta<->b ridge and lets the collapse be
# computed exactly, the way tau-slice-review.md prototyped ASIS in isolation.
# The ESS numbers are an optimistic bound on in-engine mixing; the RANKING is
# the transferable finding.
#
# Three tau chains, identical model/target, only the mean+ranef update differs:
#   CENTERED  -- beta | b (GLS on y - b_g), then b | beta, then tau | b.
#                The engine's separate f and b Gibbs blocks (chain.hpp sweep).
#   COLLAPSED -- beta | tau with b MARGINALIZED (GLS under the compound-symmetry
#                group covariance sigma^2 I + tau^2 ZZ'), then b | beta, then
#                tau | b. Draws (beta, b) as a joint block: the PRIMARY remedy.
#   ASIS      -- CENTERED plus the (tau, b) interweave (tau-slice-review.md 4c /
#                bcf-ridge-interweaving.md): eta = b/tau ancillary, redraw
#                tau | eta, y, then b <- tau eta. The SECONDARY remedy; it
#                attacks the tau-b funnel, NOT the beta-b confounding.
#
# tau is drawn by the engine's exact Makalic-Schmidt cauchy step in every arm,
# so the tau kernel is held fixed and only the confounding treatment varies.
# All three must hit the SAME posterior (a correctness self-check, printed);
# only the IACT differs. Load-independent metric; no quiet machine needed.
#
# Self-contained, seed-fixed. Run: Rscript benchmarks/R/forest-ranef-collapse-proto.R

if (!requireNamespace("coda", quietly = TRUE)) {
  stop("needs coda")
}

SIGMA <- 1.0 # known residual sd
A_CAUCHY <- 2.5 # half-Cauchy prior scale on tau (engine uses 2.5 rel.scale)
N_ITER <- 8000L
N_BURN <- 2000L

iact <- function(x) length(x) / as.numeric(coda::effectiveSize(as.numeric(x)))

## exact Makalic-Schmidt half-Cauchy tau draw (the engine's cauchy branch,
## model.hpp drawTauCauchyExactIG): xi | tau ~ IG(1, 1/tau^2 + 1/A^2);
## tau^2 | b, xi ~ IG((K+1)/2, 0.5 SS + 1/xi). IG(s, rate) = 1 / rgamma(s, rate).
draw_tau <- function(tau_cur, SS, K) {
  xi <- 1 / rgamma(1, 1, rate = 1 / tau_cur^2 + 1 / A_CAUCHY^2)
  sqrt(1 / rgamma(1, 0.5 * (K + 1), rate = 0.5 * SS + 1 / xi))
}

## Sigma^{-1} v for Sigma = sigma^2 I + tau^2 ZZ' (block compound symmetry),
## applied group-by-group via the rank-1 Woodbury downdate -- the exact form the
## in-engine collapse would need per group (design note section on the primary).
Sinv <- function(v, tau2, grp, nj) {
  gsum <- as.numeric(tapply(v, grp, sum))[grp]
  shrink <- tau2 / (SIGMA^2 * (SIGMA^2 + nj * tau2))
  v / SIGMA^2 - shrink[grp] * gsum
}

simulate <- function(n, K, s, tau_true, beta1, seed) {
  set.seed(seed)
  grp <- rep(seq_len(K), length.out = n)
  gc <- qnorm((seq_len(K) - 0.5) / K)
  x <- gc[grp] + rnorm(n, 0, s)
  b <- rnorm(K, 0, tau_true)
  y <- 0.5 + beta1 * x + b[grp] + rnorm(n, 0, SIGMA)
  list(
    y = y,
    x = x,
    grp = grp,
    K = K,
    X = cbind(1, x),
    nj = as.numeric(table(grp))
  )
}

run <- function(dat, method) {
  y <- dat$y
  X <- dat$X
  grp <- dat$grp
  K <- dat$K
  nj <- dat$nj
  n <- length(y)
  priorPrec <- diag(c(1e-4, 1e-4)) # vague beta prior
  beta <- c(0, 0)
  b <- rep(0, K)
  tau <- 0.5
  keep <- numeric(N_ITER)
  for (it in seq_len(N_ITER)) {
    if (method == "COLLAPSED") {
      ## beta | tau, marginalizing b: GLS under the group covariance.
      tau2 <- tau^2
      SiX <- apply(X, 2, Sinv, tau2 = tau2, grp = grp, nj = nj)
      Siy <- Sinv(y, tau2, grp, nj)
      prec <- crossprod(X, SiX) + priorPrec
      V <- solve(prec)
      m <- V %*% crossprod(X, Siy)
      beta <- as.numeric(m + t(chol(V)) %*% rnorm(2))
    } else {
      ## beta | b: ordinary regression of y - b_g on X (the engine's f-block).
      r <- y - b[grp]
      prec <- crossprod(X) / SIGMA^2 + priorPrec
      V <- solve(prec)
      m <- V %*% (crossprod(X, r) / SIGMA^2)
      beta <- as.numeric(m + t(chol(V)) %*% rnorm(2))
    }
    ## b | beta (conjugate group means; the engine's drawGroupEffects)
    resid <- y - as.numeric(X %*% beta)
    gsum <- as.numeric(tapply(resid, grp, sum))
    bprec <- nj / SIGMA^2 + 1 / tau^2
    bmean <- (gsum / SIGMA^2) / bprec
    b <- bmean + rnorm(K) / sqrt(bprec)
    ## tau | b (exact cauchy, held fixed across methods)
    tau <- draw_tau(tau, sum(b^2), K)
    ## ASIS (tau, b) interweave: recentre through the ancillary eta = b / tau.
    ## tau | eta, y with b = tau eta substituted is a truncated Gaussian in the
    ## ridge coefficient; here we use the exact conditional given the group means.
    if (method == "ASIS") {
      eta <- b / tau
      rbar <- gsum / nj # group means of y - X beta
      ## tau | eta ~ N(mhat, 1 / prec) truncated > 0 (half-Cauchy via the same
      ## scale mixture on the coefficient); prec = sum nj eta^2 / sigma^2 + 1/v
      v <- 1 / rgamma(1, 1, rate = (tau^2 + A_CAUCHY^2) / 2)
      prec <- sum(nj * eta^2) / SIGMA^2 + 1 / v
      mhat <- sum(nj * rbar * eta / SIGMA^2) / prec
      repeat {
        cand <- mhat + rnorm(1) / sqrt(prec)
        if (cand > 0) {
          tau <- cand
          break
        }
      }
      b <- tau * eta
    }
    keep[it] <- tau
  }
  keep[(N_BURN + 1):N_ITER]
}

## ---- confounded regime: small K, x clustered by group ---------------------
scen <- expand.grid(K = c(3L, 10L), tau_true = c(0.2), beta1 = c(3.0))
methods <- c("CENTERED", "ASIS", "COLLAPSED")
cat(sprintf(
  "n per cell = 300, s(within-group x) = 0.8, sigma = %.1f, beta1 = 3, iters = %d\n",
  SIGMA,
  N_ITER
))
cat(
  "tau kernel = exact Makalic-Schmidt cauchy in every arm; only the f/b update varies\n\n"
)
cat(sprintf(
  "%3s %8s | %-10s %-10s %-10s | %s\n",
  "K",
  "tau_true",
  "CENTERED",
  "ASIS",
  "COLLAPSED",
  "median tau (must agree)"
))
for (r in seq_len(nrow(scen))) {
  seeds <- 1:6
  ess <- setNames(numeric(length(methods)), methods)
  pm <- setNames(numeric(length(methods)), methods)
  for (m in methods) {
    ia <- numeric(length(seeds))
    md <- numeric(length(seeds))
    for (si in seq_along(seeds)) {
      dat <- simulate(
        300L,
        scen$K[r],
        0.8,
        scen$tau_true[r],
        scen$beta1[r],
        seed = 4000L + 10L * r + seeds[si]
      )
      ch <- run(dat, m)
      ia[si] <- iact(ch)
      md[si] <- median(ch)
    }
    ess[m] <- mean(ia)
    pm[m] <- mean(md)
  }
  cat(sprintf(
    "%3d %8.2f | IACT %-5.1f  IACT %-5.1f  IACT %-5.1f | %.3f / %.3f / %.3f\n",
    scen$K[r],
    scen$tau_true[r],
    ess["CENTERED"],
    ess["ASIS"],
    ess["COLLAPSED"],
    pm["CENTERED"],
    pm["ASIS"],
    pm["COLLAPSED"]
  ))
}
cat(
  "\nCOLLAPSED / CENTERED IACT ratio is the primary remedy's isolation gain;\n"
)
cat("ASIS / CENTERED is the secondary (funnel-only) gain. Median tau (heavy-\n")
cat(
  "tailed at small K, so the median not the mean) must agree across arms --\n"
)
cat("same target; a mismatch means a sampler is wrong.\n")
